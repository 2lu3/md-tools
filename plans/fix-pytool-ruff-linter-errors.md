# fix: pytool の Ruff ALL 設定による linter エラー修正

## 背景

グローバル Ruff 設定（`~/.config/ruff/ruff.toml`）で `select = ["ALL"]` が有効化されたため、`pytool/` に大量の linter エラーが出ている（初期 ~1000+）。

関連: GitHub Issue #1

## 方針

1. `pytool/pyproject.toml` に `[tool.ruff]` を追加し、リポジトリ内でも `select = ["ALL"]` を明示
2. 以下のみ ignore（修正不能 / フォーマッタ衝突 / プロジェクト方針）:
   - `COM812` / `ISC001` — formatter と衝突
   - `D203` / `D213` — 互換ルールと衝突（Ruff が自動 ignore するが明示）
   - `CPY001` — 著作権ヘッダはルート `LICENSE` で管理
   - `N813` — `import MDAnalysis as mda` は MDAnalysis の慣例
3. `extend-exclude`: `pytool/config/template`, `pytool/command/template`（生成用テンプレート）
4. `tests/**` は `S101` / `INP001` / `D` / `ANN` / `ERA001` を per-file ignore
5. 自動修正可能なものは `ruff check --fix` / `--unsafe-fixes` で適用
6. 残りはカテゴリ別に手動修正:
   - **ANN**: 型注釈追加
   - **D**: docstring 追加・整形
   - **PTH**: `os.path` / `open` → `pathlib`
   - **F401** (`__init__.py`): `__all__` で再エクスポート明示
   - **その他**: EM/TRY/FBT/S101/E501 等を個別対応
7. 既存テストが通ることを確認

## Acceptance criteria

- [x] `cd pytool && uv run ruff check pytool tests` がエラー 0
- [x] 実行可能な既存テストが通る（`test_get_box_size` / `test_reduce_dcd` / `test_rename2charmm`）
- [ ] 補足: `test_input_builder` は存在しない `pytool.input_builder` を参照（main 時点から壊れている）
- [ ] 補足: `test_log_analysis` は `tests/data/min.log` がリポジトリに無く失敗する（main 時点から）
