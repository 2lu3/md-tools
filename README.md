# md-tools

Hub repository for molecular dynamics tooling. The packages formerly
hosted here now live in their own repositories:

| Package | Repository | Install |
|---------|------------|---------|
| **dit** | [2lu3/dit](https://github.com/2lu3/dit) | `pip install git+https://github.com/2lu3/dit.git` |
| **pytool** | [2lu3/pytool](https://github.com/2lu3/pytool) | `pip install git+https://github.com/2lu3/pytool.git` |

## Migration notes

Old subdirectory installs no longer apply after this split:

```bash
# deprecated
pip install git+https://github.com/2lu3/md-tools.git#subdirectory=dit
pip install git+https://github.com/2lu3/md-tools.git#subdirectory=pytool
```

Use the repository URLs above instead.

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
