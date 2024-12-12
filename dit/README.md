# Dit

## コンセプト

* GitとDVCを統合したコマンドラインツール
* Git, DVCを直接扱うことも可能
* 特定のディレクトリのみを対象にDVCを使うことで、ストレージを節約することができる

## コマンド

### `add <path>`

* `dvc add` 実行後に`git add` を実行する
* scopeに自動追加される

### `commit -m "message"`

* `git commit` を実行する

### `push`

* `dvc push` 実行後に `git push` を実行する

### `scope`

* scopeを追加/削除する
* インタラクティブに操作する

#### `scope add <path>`

* 指定したディレクトリをscopeに追加する

#### `scope remove <path>`

* 指定したディレクトリをscopeから削除する


#### `scope list`

* scopeに含まれるディレクトリを表示する


