# マウスTtn遺伝子のスプライシングに関する情報の可視化
## 下記の情報についての可視化を実装
1. 遺伝子上のTCTT配列の分布を可視化
2. 異なる長さのイントロンの分布を可視化

## visualization.ipynbで、可視化された情報を確認
NCBIののマウスTtn遺伝子のgenbankファイルよりデータを読み込み可視化を行う

## visualization.ipynbの内容
1. マウスTtnの異なるバリアントとのエクソン、イントロン領域におけるTCTT配列の分布を確認
2. エクソン、イントロンに分け、各イントロンの塩基あたりのTCTT出現率を比較
3. 異なる生物種のバリアントごとに見られイントロン領域のパターンを可視化
4. エクソン数の大きい遺伝子の各バリアントのイントロン領域のパターンを可視化

## 1. 遺伝子上のTCTT配列の分布の可視化に関するコマンド
#### 準備
```
gbk = Seq_count() # クラスのインスタンス化
gbk.read_gbk('data/gbk/mouse_ttn.gb')　# gbkファイルの読み込み
gbk.tutorial() # tutorialの表示
```

#### スプライスバリアント関連
```
gbk.get_mrna_ids() # バリアントの表示
gbk.set_mrna_id('NM_001385708.1') # トランスクリプトバリアントの設定
```

#### TCTT配列の設定および分布の可視化
```
gbk.set_interest_seq('TCTT') # TCTT配列の設定
gbk.heatmap_hist() # 着目する配列の分布の可視化
gbk.save_fig()　# ヒートマップをfig/に保存
```

#### イントロン領域のTCTT配列の分布の棒グラフ
```
gbk.intron_bar() # 各イントロン内のTCTT配列の棒グラフ
gbk.intron_bar_base() # 各イントロン内の100塩基あたりのTCTT配列の棒グラフ
gbk.interest_seq_count_edge() #  各イントロンの両端の５０塩基に存在するTCTT配列の棒グラフ
```

#### イントロン領域のTCTT配列の分布の棒グラフ
```
gbk.exon_bar() # 各エクソン内のTCTT配列の棒グラフ
gbk.exon_bar_base() # 各エクソン内の100塩基あたりのTCTT配列の棒グラフ
```

## 2. 異なる長さのイントロンの分布の可視化に関するコマンド
#### NCBI情報
RNA-seq intron features, aggregate (filterd)のファイルを取得<br>
[**mouse Ttn**](https://www.ncbi.nlm.nih.gov/nuccore/NC_000068.8?report=graph&from=76492536&to=76854687&strand=true&app_context=Gene&assm_context=GCF_000001635.27)<br>
[**rat Ttn**](https://www.ncbi.nlm.nih.gov/nuccore/NC_051338.1?report=graph&from=61611559&to=61965783&strand=true&app_context=Gene&assm_context=GCF_015227675.2)<br>
[**rabbit Ttn**](https://www.ncbi.nlm.nih.gov/nuccore/NC_067380.1?report=graph&from=55045025&to=55398485&app_context=Gene&assm_context=GCF_009806435.1)<br>
[**human TTN**](https://www.ncbi.nlm.nih.gov/nuccore/NC_000002.12?report=graph&from=178483773&to=178849637&strand=true&app_context=Gene&assm_context=GCF_000001405.40)<br>
[**chimp TTN**](https://www.ncbi.nlm.nih.gov/nuccore/NC_072400.1?report=graph&from=80513939&to=80881842&strand=true&app_context=Gene&assm_context=GCF_028858775.1)<br>
[**gorilla TTN**](https://www.ncbi.nlm.nih.gov/nuccore/NC_073226.1?report=graph&from=69432214&to=69798206&strand=true&app_context=Gene&assm_context=GCF_029281585.1)
#### 準備
```
pile = PileUp() # クラスのインスタンス化
pile.set_csv(f"data/intron/mouse_ttn_intron.tsv") # NCBIから取得したイントロンに関するデータを読み込み
```
#### 異なる長さのイントロンの分布を表示
```
pile.set_title(f"Mouse Ttn intron gap distribution")　# タイトルの設定
pile.show(height=300) # グラフの表示
```




## 種間の保存性の可視化
[**UCSC genome browser**](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A178525989%2D178807423&hgsid=1735195872_llHo1dKi3SVLwkPP3OCN3pMDqYzb)