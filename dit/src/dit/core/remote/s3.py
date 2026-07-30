from __future__ import annotations

from pathlib import Path

import boto3
from botocore.config import Config as BotoConfig
from botocore.exceptions import ClientError

from dit.core.config import RemoteConfig, S3_REQUEST_TIMEOUT_SEC
from dit.core.hasher import strip_hash_prefix
from dit.core.remote.base import Remote

MULTIPART_THRESHOLD_BYTES = 64 * 1024 * 1024
MULTIPART_CHUNKSIZE_BYTES = 16 * 1024 * 1024
MAX_CONCURRENCY = 4


class S3Remote(Remote):
    def __init__(self, remote: RemoteConfig) -> None:
        self.remote = remote
        self._client = boto3.client(
            "s3",
            endpoint_url=remote.endpoint_url,
            config=BotoConfig(
                connect_timeout=S3_REQUEST_TIMEOUT_SEC,
                read_timeout=S3_REQUEST_TIMEOUT_SEC,
                retries={"max_attempts": 3, "mode": "standard"},
            ),
        )
        self._transfer = boto3.s3.transfer.TransferConfig(
            multipart_threshold=MULTIPART_THRESHOLD_BYTES,
            multipart_chunksize=MULTIPART_CHUNKSIZE_BYTES,
            max_concurrency=MAX_CONCURRENCY,
        )

    def object_key(self, content_hash: str) -> str:
        digest = strip_hash_prefix(content_hash)
        prefix = self.remote.prefix.rstrip("/")
        base = f"files/blake3/{digest[:2]}/{digest[2:]}"
        return f"{prefix}/{base}" if prefix else base

    def exists(self, content_hash: str) -> bool:
        key = self.object_key(content_hash)
        try:
            self._client.head_object(Bucket=self.remote.bucket, Key=key)
            return True
        except ClientError as exc:
            code = exc.response.get("Error", {}).get("Code")
            if code in {"404", "NoSuchKey", "NotFound"}:
                return False
            raise RuntimeError(
                f"failed to head s3://{self.remote.bucket}/{key}: {exc}"
            ) from exc

    def upload(self, local_path: Path, content_hash: str) -> None:
        key = self.object_key(content_hash)
        if self.exists(content_hash):
            return
        try:
            self._client.upload_file(
                Filename=str(local_path),
                Bucket=self.remote.bucket,
                Key=key,
                Config=self._transfer,
            )
        except ClientError as exc:
            raise RuntimeError(
                f"failed to upload {local_path} to s3://{self.remote.bucket}/{key}: {exc}"
            ) from exc

    def download(self, content_hash: str, local_path: Path) -> None:
        key = self.object_key(content_hash)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = local_path.with_suffix(local_path.suffix + ".ditdownload")
        try:
            self._client.download_file(
                Bucket=self.remote.bucket,
                Key=key,
                Filename=str(tmp),
                Config=self._transfer,
            )
            tmp.replace(local_path)
        except ClientError as exc:
            if tmp.exists():
                tmp.unlink(missing_ok=True)
            raise RuntimeError(
                f"failed to download s3://{self.remote.bucket}/{key} to {local_path}: {exc}"
            ) from exc

    def delete(self, content_hash: str) -> None:
        key = self.object_key(content_hash)
        try:
            self._client.delete_object(Bucket=self.remote.bucket, Key=key)
        except ClientError as exc:
            raise RuntimeError(
                f"failed to delete s3://{self.remote.bucket}/{key}: {exc}"
            ) from exc

    def list_hashes(self) -> list[str]:
        prefix = self.remote.prefix.rstrip("/")
        list_prefix = f"{prefix}/files/blake3/" if prefix else "files/blake3/"
        hashes: list[str] = []
        continuation: str | None = None
        while True:
            kwargs = {"Bucket": self.remote.bucket, "Prefix": list_prefix}
            if continuation:
                kwargs["ContinuationToken"] = continuation
            try:
                resp = self._client.list_objects_v2(**kwargs)
            except ClientError as exc:
                raise RuntimeError(
                    f"failed to list s3://{self.remote.bucket}/{list_prefix}: {exc}"
                ) from exc
            for obj in resp.get("Contents") or []:
                key = obj["Key"]
                digest = _digest_from_key(key, list_prefix)
                if digest:
                    hashes.append(f"blake3:{digest}")
            if not resp.get("IsTruncated"):
                break
            continuation = resp.get("NextContinuationToken")
        return hashes


def _digest_from_key(key: str, list_prefix: str) -> str | None:
    if not key.startswith(list_prefix):
        return None
    rest = key[len(list_prefix) :]
    parts = rest.split("/")
    if len(parts) != 2:
        return None
    return parts[0] + parts[1]


def open_remote(remote: RemoteConfig) -> S3Remote:
    return S3Remote(remote)
