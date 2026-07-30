---
name: dit standalone MD VCS
overview: 既存の dit を DVC 非依存のスタンドアロンな大容量ファイル管理ツールに作り直す。管理対象は dit.toml で宣言し、pre-commit フックが自動でポインタを更新する。日常操作は dit sync（①.dcd/.dit 一致確認と新しい方の採用 ②scope に合わせてローカル実体を置く/消す）に集約する。孤児リモート削除は --prune-remote 明示時のみ、自動 fetch 後に実行する。
todos:
  - id: issue-and-plan
    content: GitHub issue を作成し、plans/feat-dit-standalone.md をコミット。feature ブランチを作成（origin/main を fetch して最新から分岐）
    status: pending
  - id: pyproject
    content: dit/pyproject.toml を poetry から uv + hatchling(PEP 621) に移行し、src layout (packages = ["src/dit"]) を設定。click/boto3/blake3/pathspec/loguru/natsort/alive-progress を uv add。dev に pytest/moto を追加
    status: pending
  - id: src-layout
    content: 既存の dit/dit/ を dit/src/dit/ に git mv し、パッケージが src layout でインポートできることを確認
    status: pending
  - id: core-repo-config
    content: core/repo.py（git ルート・dit.toml・.dit 解決）と core/config.py（dit.toml の [remote]/[track] 読み込み、認証情報は環境変数から）を実装
    status: pending
  - id: core-tracker
    content: core/tracker.py を実装。pathspec の gitwildmatch で dit.toml のパターンを解釈し、.git/ と .gitignore を尊重してワークツリーを走査。パターン解釈のテストを追加
    status: pending
  - id: core-hash-index
    content: core/hasher.py（BLAKE3 並列ハッシュ）と core/index.py（sqlite の size/mtime_ns/inode インデックス、stat 一致時はハッシュ省略）を実装しテストを追加
    status: pending
  - id: core-pointer-scope
    content: core/pointer.py（*.dit の TOML read/write と .gitignore 追記）と core/scope.py（.dit/scope.toml）を実装
    status: pending
  - id: core-remote-s3
    content: core/remote/s3.py を実装。content-addressed キー、head_object による存在確認、TransferConfig によるマルチパート並列転送。moto でテスト
    status: pending
  - id: cmd-basic
    content: init（既定パターン入り dit.toml の生成 + フック設置）/ add（引数なし。リポジトリ全体を走査してポインタ更新、--quiet と --prune のみ）/ status を実装
    status: pending
  - id: githook
    content: core/githook.py と hooks/pre-commit テンプレート、dit hook install|uninstall|status を実装。既存フックの上書き回避、core.hooksPath 対応、dit 未インストール時の明示エラー、部分コミット時の警告を含める
    status: pending
  - id: cmd-transfer
    content: push / pull（低レベル）/ sync（方針1+2、--dry-run、--prune-remote 時は自動 git fetch --all --prune 後に全 ref 走査で孤児削除）/ scope コマンドを実装
    status: pending
  - id: cleanup-docs
    content: 旧 DVC 依存コード（model/dvc.py, gui/gui.py, model/command.py, model/extension.py, command/extension.py）を削除し、dit/README.md とリポジトリ README.md を dit.toml / sync ベースの新仕様に更新
    status: pending
isProject: true
---

# dit を DVC 非依存の MD 向けバージョン管理ツールに作り直す

## 決定事項

- ベース: 既存の [dit/](dit/) を作り直す。DVC には一切依存しない
- パッケージ: **src layout**（`dit/src/dit/`）
- リモート: S3 互換（MinIO）のみを v1 で実装。boto3 を使用
- キャッシュ: **ローカルキャッシュを持たない**。ワークツリーのファイルが実体、リモートが真実
- 変更検知: `.dit/index.db`（sqlite）に `size / mtime_ns / inode` を記録し、一致すればハッシュを再計算しない。変わった時だけ BLAKE3 で全読み
- Git に記録するもの: ファイル単位ポインタ `foo.dcd` → `foo.dcd.dit`
- **管理対象は宣言的に決める**: リポジトリルートの `dit.toml`（`.gitignore` 書式、拡張子・ディレクトリ・`!` 除外）
- **git との連携は pre-commit フック**: `git commit` のたびに `dit add` を自動実行。`add` に path 引数は無い
- v1 コマンド: `init` / `add` / `status` / `push` / `pull` / `sync` / `scope` / `hook`
- **`gc` / `fsck` は置かない**（`sync` と `sync --dry-run` が担う）

## `sync` の方針（2つだけ）

1. **ローカルとリモートの新しい方を採用する**（必要なら push / pull）。`.dit` もそれに合わせて正しい情報にする。**.dcd と .dit が両方ある場合は必ず一致確認**する（stat で変わっていなければ index のハッシュで足りる。変わっていれば BLAKE3 で再計算）。不一致なら mtime で新しい方を採用
2. **scope 内ならローカルに実体を置き、scope 外ならローカルの実体を消す**（消す前にリモート存在を確認。無ければ削除せず error）

### 状態処理（方針から導出）

| .dit | .dcd | remote | 動作 |
|------|------|--------|------|
| あり | あり（ハッシュ一致） | あり | 一致確認 OK。scope 内 → 正常 / scope 外 → ローカル削除（方針2） |
| あり | あり（ハッシュ不一致） | — | 一致確認で不一致 → mtime で新しい方を採用し `.dit` を更新。必要なら push / pull。その後方針2 |
| あり | あり | なし | push（方針1）。その後方針2（scope 外なら削除） |
| あり | なし | あり | scope 内 → pull / scope 外 → 正常（方針2） |
| あり | なし | なし | **error**: ファイルが見つかりません |
| なし | あり | あり | **warning**: `hoge.dcd.dit` がありません |
| なし | あり | なし | **warning**: untracked |
| なし | なし | あり | 既定では触らない。`--prune-remote` 時のみ削除（下記） |
| なし | なし | なし | 対象外（観測できない） |

### 孤児リモート削除（`--prune-remote`）

ローカル ref だけでは fetch していない remote ブランチの参照が見えない。

- **既定では孤児削除しない**
- `dit sync --prune-remote` の明示時だけ実行する
- 実行時は **自動で `git fetch --all --prune` してから**、全ローカル ref（fetch 後）のポインタを走査し、未参照オブジェクトだけ削除する
- `--dry-run` と併用可

## `dit.toml`（リポジトリルート、git 管理）

```toml
[remote]
url = "s3://my-bucket/md-project"
endpoint_url = "https://minio.example.com"

[track]
patterns = [
  "*.dcd",
  "*.dvl",
  "*.rst",
  "*.npy",
  "*.pkl",
  "*.tar",
  "data/**/out/",
  "!data/00_scratch/",
]
```

パターン解釈は pathspec の `gitwildmatch`。認証情報は環境変数のみ（`AWS_ACCESS_KEY_ID` 等）。

## ディレクトリ構成

```
dit/
├── pyproject.toml               # uv + hatchling, src layout
├── src/
│   └── dit/
│       ├── main.py
│       ├── core/
│       │   ├── repo.py
│       │   ├── config.py
│       │   ├── tracker.py
│       │   ├── pointer.py
│       │   ├── hasher.py
│       │   ├── index.py
│       │   ├── scope.py
│       │   ├── remote/s3.py
│       │   └── githook.py
│       ├── hooks/pre-commit
│       └── command/             # init/add/status/push/pull/sync/scope/hook
└── tests/test_*.py
```

### `.dit/`（git 管理外、マシンローカル）

- `index.db`: `path, size, mtime_ns, inode, hash, pushed_at`
- `scope.toml`: このマシンで実体を持ちたいディレクトリ
- `.gitignore`: `*`

### リモートキー

`<prefix>/files/blake3/<hash[:2]>/<hash[2:]>`

## 各コマンド

- `init` — 既定パターン入り `dit.toml` + `.dit/` + pre-commit フック
- `add` — **引数なし**。`dit.toml` に一致するファイルを走査してポインタ更新。`--quiet` / `--prune` のみ。通常はフックから呼ばれる
- `hook install|uninstall|status`
- `status` — `M` / `?` / `↓` / `↑` 等
- `push` / `pull` — sync から呼ばれる低レベル操作
- `sync` — 日常の主入口。方針1+2。`--dry-run` / `--prune-remote`
- `scope add|remove|list` — このマシンで実体を持ちたいディレクトリ

## scope と dit.toml

| | 決めること | 保存先 | 共有 |
|---|---|---|---|
| `dit.toml` | 何を dit の管理対象にするか | リポジトリルート | git で共有 |
| scope | このマシンのディスクに実体を置きたいのはどこか | `.dit/scope.toml` | 共有しない |

## v1 で扱わないこと

- 既存 `.dvc` からの移行（`dit migrate`）
- SSH/rsync リモート、チャンク分割による差分転送
- パイプライン実行（`dvc repro` 相当）、pytool の step レイアウト連携
- `gc` / `fsck` コマンド

## 進め方

GitHub issue → `plans/feat-dit-standalone.md` をコミット → feature ブランチで実装 → PR。
core（hash / index / pointer / remote）を先に固め、テストを付けてからコマンド層を載せる。
