import plotly.express as px
import pandas as pd
import numpy as np

class PileUp():

    def set_csv(self, csv):
        self.df = pd.read_csv(csv, sep='\t', header=None, usecols=[1,2])
    
    def set_title(self, title):
        self.set_title = title

    def show(self, height=600):
        self.height = height
        numlst = [[i,j] for i,j in zip(self.df[1], self.df[2])]
        # それぞれのリストの長さを追加
        numlst = [[i,j,j-i] for i,j in numlst]

        # 始点の位置とその長さ順にソート
        numlst.sort(key=lambda x :(x[0], -x[2]))

        # numlstの表示
        # [始点, 終点, 長さ]
        #print(numlst)

        # ソート後のnumlstから長さの情報を削除
        numlst = [[i,j] for i,j,k in numlst]
        #print(numlst)

        # result辞書を作成
        # keyに階層情報
        # valueに始点、終点情報を代入
        result = {0:[numlst[0]]}

        # result辞書の更新
        # 一番下の階層の最終リストの終点が、入力した始点の値より小さい場合は、同一の階層のvalueにappendする
        # そうでない場合は、一つ上の階層で同様の処理をする
        # その処理が行えない場合は、新たな階層のkeyを作成し、そこに始点と終点のリストをappendする
        for i in numlst[1:]:
            for j in range(len(result.keys())):
                if result[j][-1][1] < i[0]:
                    result[j].append(i)
                    break
            else:
                result.setdefault(max(result.keys())+1,[i])           
        #print(result)

        # figureの作図
        # [始点,終点,None, 始点, 終点, None, .....]
        # [階層#, 階層#, ......]
        # 上記の一次元のリストを作成し、それぞれx,yの値としてpx.lineで作図
        x_arr = np.array([])
        y_arr = np.array([])
        for i in range(len(result.keys())):
            x_arr = np.concatenate([x_arr, np.array( [i + [None] for i in result[i]]).flatten()])
            y_arr = np.concatenate([y_arr, np.full(3*len(result[i]), i)])
        fig = px.line(x=x_arr, y=y_arr)
        fig.update_traces(line=dict(width=10))
        fig.update_layout(title=self.set_title, title_x=0.5, width=1000, height=self.height)
        fig.show()