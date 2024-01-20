import openpyxl as px

#ファイルを開く
# path = r"C:\Users\watanabe_ryunosuke\PythonPrograms\Medical_Stem_Modeling\Sub\ds.dat"
path = "Input\ds\ds.dat"
file = open(path)

#すべての行をリストとして読み込み
lines = file.readlines()
#print(lines)

#ファイルを閉じる
file.close()

#ブックを新規作成。Workbookの'w'は必ず大文字。
wb = px.Workbook()

#出力ファイル名を指定
filename = "Output/preprocess.csv"

#アクティブなシートを取得
ws = wb.active


"""
以下、テキストファイルから取得した値を、想定したデザイン通りになるように
セルに書き込んでいく。
"""
#listを取得するための変数
l = 0

#セルへの書き込み
for i in range(1,1772+1):
        blank = lines[l].split()
        for j in range(0, len(blank) ):
            ws.cell(row=i, column=j+1).value = float(blank[j])
        l += 1
        print(i)

#ブックを保存する
wb.save(filename)
