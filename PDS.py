from modules import CrossingNumberAlgorithm
from modules import CheckDistance
from scipy.interpolate import Rbf
import numpy as np
import csv
import time
import random

# define
MAXIMUM_NUMBER_OF_SEARCHES = 1000 # 点が連続でN回生成できなかったら終了
MAXIMUM_NUMBER_OF_POINTS = 1000 # 物体内部最大生成点数
PITCH_SURFACE = 1 # 物体表面に点群を生成するときに用いる．
PITCH_MABIKI = 3 # 物体表面に生成した点から間引きを行う際のPDSの最小点間距離
ALLOWABLE_STRESS = 186 #チタン合金．降伏強さ930MPa．安全率5
RADIUS = 0.5 # ラティス半径[mm]

# Inputファイル
input_path = 'Input/Stress/HourGrass2_stress.csv' # ANSYSのデータファイル
mesh_data = 'Input/Mesh_Data/Hourgrass2_mesh.stl' # 物体の表面形状データ。

# Outputファイル
surface_Points_path = 'Output/PDS/surface6.ply' # 表面点のみ表示
inner_Points_path = 'Output/PDS/inner4.ply' # 内部点のみ表示
result_ply_path = 'Output/PDS/pds4.ply'
result_csv_path = 'Output/PDS/pds4.csv'

def main():
    surface_Points = [] # 表面確定点格納用
    inner_Points = [] # 内部確定点格納用

    # ANSYS上の点群を取得し座標値を取得
    points = np.loadtxt(input_path, delimiter=',')
    print('points = ',points[0:3])
        
    rbf = Rbf(points[:,1], points[:,2], points[:,3], points[:,4], function='multiquadric')
    
    # PDSでの点の生成範囲の設定
    x_max = np.amax(points[:,1])
    x_min = np.amin(points[:,1])
    y_max = np.amax(points[:,2])
    y_min = np.amin(points[:,2])
    z_max = np.amax(points[:,3])
    z_min = np.amin(points[:,3])
    
    # PDSでの点の生成範囲の表示
    print("x_max = ", x_max, "x_min = ", x_min)
    print("y_max = ", y_max, "y_min = ", y_min)
    print("z_max = ", z_max, "z_min = ", z_min)

    #time.sleep(10)

    # 交差数判定法
    CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm(mesh_data)

    # PDS用
    CD = CheckDistance.CheckDistance(ALLOWABLE_STRESS, RADIUS)

    #物体表面上でPDS
    CNA.surface_generate(surface_Points, PITCH_SURFACE, PITCH_MABIKI)
    #重複した座標を削除
    surface_Points, _ = np.unique(surface_Points, return_index=True, axis=0)


    num = 0
    while num < MAXIMUM_NUMBER_OF_SEARCHES:
        if len(inner_Points) >= MAXIMUM_NUMBER_OF_POINTS:
            break

        flg_P = False
        while not flg_P: # 物体内部の点を生成するまでループ。内部であればflg_PはTrueとなる。
            # ランダムな点を生成
            pds_x = random.uniform(x_min, x_max)
            pds_y = random.uniform(y_min, y_max)
            pds_z = random.uniform(z_min, z_max)
            candidate_point = [pds_x, pds_y, pds_z]
    
            # 物体内部の点か判定
            flg_P = CNA.cramer(candidate_point)

        new_stress = rbf(candidate_point[0], candidate_point[1], candidate_point[2])
        # 点間距離内に他の点が含まれているか否かを判定
        flg = CD.check_distance(inner_Points, candidate_point, new_stress)

        # 点間距離内に他の点が存在しないとき候補点を確定点に追加
        if flg:
            inner_Points.append(candidate_point)
            num = 0
            print('num : ', num, 'inner_Points : ', len(inner_Points))
        # 点間距離内に他の点が存在するとき
        else :
            num = num + 1
            print('num : ', num)


    # ply にPDSの表面点のみ結果出力
    inner_Points, _ = np.unique(inner_Points, return_index=True, axis=0) # 重複した座標を削除
    with open(surface_Points_path, 'w', newline="") as f:
        f.write('ply\n')
        f.write('format ascii 1.0\n')
        f.write('element vertex '+str(len(surface_Points))+'\n')
        f.write('property double x\n')
        f.write('property double y\n')
        f.write('property double z\n')
        f.write('end_header\n')
        for ele in surface_Points:
            f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')

    # ply にPDSの内部点のみ結果出力
    inner_Points, _ = np.unique(inner_Points, return_index=True, axis=0) # 重複した座標を削除
    with open(inner_Points_path, 'w', newline="") as f:
        f.write('ply\n')
        f.write('format ascii 1.0\n')
        f.write('element vertex '+str(len(inner_Points))+'\n')
        f.write('property double x\n')
        f.write('property double y\n')
        f.write('property double z\n')
        f.write('end_header\n')
        for ele in inner_Points:
            f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')

    # ply にPDSの結果出力
    fixed_Points = np.concatenate([surface_Points, inner_Points]) # 結合
    fixed_Points, _ = np.unique(fixed_Points, return_index=True, axis=0) # 重複した座標を削除
    with open(result_ply_path, 'w', newline="") as f:
        f.write('ply\n')
        f.write('format ascii 1.0\n')
        f.write('element vertex '+str(len(fixed_Points))+'\n')
        f.write('property double x\n')
        f.write('property double y\n')
        f.write('property double z\n')
        f.write('end_header\n')
        for ele in fixed_Points:
            f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')

    # csv にPDSの結果出力
    with open(result_csv_path, 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(fixed_Points)

    # 点の個数を出力
    print("fixed_Points:   ", len(fixed_Points))
    print("surface_Points:   ", len(surface_Points))
    print("inner_Points:    ", len(inner_Points))
    


if __name__ == '__main__':
    start = time.time()  # 時間計測用
    main() # 点群生成
    
    # 時間計測結果の表示
    elapsed_time = time.time() - start 
    print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")

