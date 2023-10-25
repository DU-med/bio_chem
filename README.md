# 生化学講座プロジェクト

## 環境設定

### Mambaのインストール
下記のリポジトリから、mambaのインストーラを取得<br>
https://github.com/conda-forge/miniforge#mambaforge

**例：ubuntu(linux, Intel CPU)にインストールする場合**<br>
```
wget -c https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
```

### 仮想環境の構築
bioという名前の仮想環境を作成し、pythonのバージョンを3.10にする
```
mamba create -n bio python=3.10
```