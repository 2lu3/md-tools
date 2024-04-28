# md-tools

This repository contains some useful tools to performe molecular dynamics simulation.


## pytool

```
pip install git+https://github.com/2lu3/md-tools.git#subdirectory=pytool
```

## DVC

```
pipx install dvc

dvc init

dvc remote add -d minio s3://name
dvc remote modify minio endpointurl URL
```

```
dvc add --glob "**/*.dcd"
dvc add --glob "**/*.dvl"
dvc add --glob "**/*.rst"
```

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
