# md-tools

This repository contains some useful tools to perform molecular dynamics simulation.

## pytool

```
pip install git+https://github.com/2lu3/md-tools.git#subdirectory=pytool
pip install git+https://github.com/2lu3/md-tools.git#subdirectory=dit
```

## dit

MD 向けの大容量ファイル版管理（DVC 非依存）。詳細は [dit/README.md](dit/README.md)。

```bash
pip install git+https://github.com/2lu3/md-tools.git#subdirectory=dit

dit init --remote-url s3://bucket/prefix --endpoint-url https://minio.example.com
# edit dit.toml [track].patterns as needed
dit scope add data/07_production
dit sync
```

`git commit` 時に pre-commit フックが `dit add` を実行し、`*.dit` ポインタを自動ステージする。

## .env and .envrc

.env

```
AWS_ACCESS_KEY_ID=""
AWS_SECRET_ACCESS_KEY=""
FUGAKU_USER_ID=""
```

.envrc

```
dotenv
export PATH="$PATH:$(pwd)/software/genesis/bin"
export PATH="$PATH:$(pwd)/software/genesis-2.1.0-cpu/bin"
export PATH="$PATH:$(pwd)/tool/bin"
```
