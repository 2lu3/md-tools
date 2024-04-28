#!/bin/bash

# benchmarksディレクトリ以下を検索対象とする
search_dir="benchmarks"

# 指定されたパターンのファイルを検索
files=$(find $search_dir -type f \( -path "*/out/*.dcd" -o -name "*/out/*.rst" \))

if [ -z "$files" ]; then
  echo "No files found to delete."
  exit 0
fi

# 検索結果を表示
echo "The following files will be deleted:"
echo "$files"

# ユーザーに削除の確認を取る
echo "Do you want to delete these files? [y/n]"
read answer

if [ "$answer" = "y" ]; then
  # ファイル削除
  echo "$files" | xargs rm
  echo "Files deleted."
else
  echo "File deletion cancelled."
fi

