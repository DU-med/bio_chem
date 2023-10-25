# マウスTtn遺伝子上のTCTT配列の分布を可視化
## visualization.ipynbで、可視化された情報を確認
NCBIののマウスTtn遺伝子のgenbankファイルよりデータを読み込み可視化を行う

## visualization.ipynbの内容
1. マウスTtnの異なるバリアントとのエクソン、イントロン領域におけるTCTT配列の分布を確認
2. エクソン、イントロンに分け、各イントロンの塩基あたりのTCTT出現率を比較
3. 異なる生物種のバリアントごとに見られイントロン領域のパターンを可視化
4. エクソン数の大きい遺伝子の各バリアントのイントロン領域のパターンを可視化

## 基本的なコマンドの説明
#### クラスのインスタンス化
```
gbk = Seq_count()
```
#### gbkファイルの読み込み
```
gbk.read_gbk('data/gbk/mouse_ttn.gb')
```
#### tutorialの表示
```
gbk.tutorial()
```
#### バリアントの表示
```
gbk.get_mrna_ids()
```
#### トランスクリプトバリアントの設定
```
gbk.set_mrna_id('NM_001385708.1')
```
#### 着目する配列の設定
```
gbk.set_interest_seq('TCTT')
```
#### TCTT配列の分布の可視化
```
gbk.heatmap_hist()
```
#### figureの保存
```
gbk.save_fig()
```