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
mamba create -n bio python=3.10 -y
```
仮想環境の起動
```
mamba activate bio
```

### 基本的なライブラリのインストール
```
mamba install -c anaconda jupyter -y
mamba install -c anaconda pandas -y
mamba install -c anaconda seaborn -y
mamba install -c anaconda scipy -y
mamba install -c anaconda scikit-learn -y
mamba install -c anaconda fastcluster -y
mamba install -c anaconda statsmodels -y
mamba install -c conda-forge biopython -y
mamba install -c plotly plotly -y
mamba install -c bioconda logomaker -y

```

### リポジトリのクローン
```
git clone https://github.com/utsumidaisuke/bio_chem.git
```

### jupyter notebookの起動
```
jupyter notebook
```
