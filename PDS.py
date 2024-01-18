# 目的：ドロネー分割を節点番号表記に変更．その後，PDSによる点群生成．
# Input：ドロネー分割の結果（座標値のみ：③の結果）＆　②の結果
# Output：PDSの点群座標.plyファイル，PDSの点群座標.csvファイル
# ver2との違い：trimeshをCNAで用いるため，SPLITを消した．またCOEFFICIENT_OF_LONGも消した

from modules import Point
from modules import CrossingNumberAlgorithm
from modules import CheckDistance
from scipy.interpolate import Rbf
import numpy as np
import csv
import time
import random

# define
MAXIMUM_NUMBER_OF_SEARCHES = 800 # 点が連続でN回生成できなかったら終了
MAXIMUM_NUMBER_OF_POINTS = 600 # 物体内部最大生成点数
PITCH = 1 # kikalabの物体表面に点群を生成するときに用いる．
#RATE_OF_THINNINGS = 0.05 # 間引き後の点の割合．例：0.05→5%
ALLOWABLE_STRESS = 186 #チタン合金．降伏強さ930MPa．安全率5
PDS_PITCH = 3 # 物体表面に生成した点から間引きを行う際のPDSの最小点間距離



# Inputファイル
# input_path = './Input/Column10_0615.csv' # ANSYSのデータファイル
input_path = './Input/cube_50x50mm.csv' # ANSYSのデータファイル
mesh_data = 'Input/Mesh_Data/cube_50x50mm_mesh.stl' # 物体の表面形状データ。

# Outputファイル
resultIn_ply_path = 'Output/result_main/resultIn.ply' # 内部点のみ表示
result_ply_path = 'Output/result_main/result.ply'
result_csv_path = 'Output/result_main/result.csv'
biharmonic_path = 'Output/Biharmonic/biharmonic_param.csv'


def main():
    points_obj_list = []  # Pointオブジェクトを保持
    fixed_points = [] # 確定点格納用

    # ANSYS上の点群を取得し座標値を取得
    points = np.loadtxt(input_path, delimiter=',')
    print('points = ',points[0:3])
    
    # FEMの節点の読み込み
    print('FEM節点を読み込んでいます')
    for i in range(len(points)):
        point = Point.Point(points[i]) # クラスに格納
        point.system_guid_obj_to_coordinate()
        points_obj_list.append(point)
    
    
    rbf = Rbf(points[:,1], points[:,2], points[:,3], points[:,4], function='multiquadric')
    print('rbf:   ', rbf(points[10,1], points[10,2], points[10,3]))
    print('stress:   ', points[10,4])
    # PDSでの点の生成範囲の設定
    x_max = points[0,1]
    x_min = points[0,1]
    y_max = points[0,2]
    y_min = points[0,2]
    z_max = points[0,3]
    z_min = points[0,3]
    
    for point in points:
        if point[1] > x_max:
            x_max = point[1]
        if point[1] < x_min:
            x_min = point[1]
        
        if point[2] > y_max:
            y_max = point[2]
        if point[2] < y_min:
            y_min = point[2]

        if point[3] > z_max:
            z_max = point[3]
        if point[3] < z_min:
            z_min = point[3]
    
    # PDSでの点の生成範囲の表示
    print("x_max = ", x_max, "x_min = ", x_min)
    print("y_max = ", y_max, "y_min = ", y_min)
    print("z_max = ", z_max, "z_min = ", z_min)
    # print('super_box : ', super_box)

    # 交差数判定法
    CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm(mesh_data)


    # PDS用
    CD = CheckDistance.CheckDistance(ALLOWABLE_STRESS)

    num = 0
    while num < MAXIMUM_NUMBER_OF_SEARCHES:
        if len(fixed_points) >= MAXIMUM_NUMBER_OF_POINTS:
            break

        flg_P = False
        while not flg_P: # 物体内部の点を生成するまでループ。内部であればflg_PはTrueとなる。
            # ランダムな点を生成
            pds_x = random.uniform(x_min, x_max)
            pds_y = random.uniform(y_min, y_max)
            pds_z = random.uniform(z_min, z_max)
            pds_point = [pds_x, pds_y, pds_z]
            #print("pds_point : ",pds_point)
    
            # 物体内部の点か判定
            flg_P = CNA.cramer(pds_point)

        # 生成点の座標情報をPointに格納し候補点にする
        candidate_point = Point.Point(pds_point)
        candidate_point.pds_coordinate()
        new_stress = rbf(candidate_point.x, candidate_point.y, candidate_point.z)
        # 点間距離内に他の点が含まれているか否かを判定
        flg = CD.check_distance(fixed_points, candidate_point, new_stress)

        # 点間距離内に他の点が存在しないとき候補点を確定点に追加
        if flg:
            fixed_points.append(pds_point)
            num = 0
            print('num : ', num, 'fixed_points : ', len(fixed_points))
        # 点間距離内に他の点が存在するとき
        else :
            num = num + 1
            print('num : ', num)


    # ply にPDSの内部点のみ結果出力
    #重複した座標を削除
    fixed_Inpoints = fixed_points
    fixed_Inpoints, _ = np.unique(fixed_points, return_index=True, axis=0)
    with open(resultIn_ply_path, 'w', newline="") as f:
        f.write('ply\n')
        f.write('format ascii 1.0\n')
        f.write('element vertex '+str(len(fixed_Inpoints))+'\n')
        f.write('property double x\n')
        f.write('property double y\n')
        f.write('property double z\n')
        f.write('end_header\n')
        for ele in fixed_Inpoints:
            f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')

    #max，minのdensityを表示
    CD.print_max_min_density()

    #物体表面上でPDS
    CNA.surface_kikalab(fixed_points, PITCH, PDS_PITCH)

    #重複した座標を削除
    fixed_points, _ = np.unique(fixed_points, return_index=True, axis=0)

    # ply にPDSの結果出力
    print("fixed_points  = ", len(fixed_points), "個")
    with open(result_ply_path, 'w', newline="") as f:
        f.write('ply\n')
        f.write('format ascii 1.0\n')
        f.write('element vertex '+str(len(fixed_points))+'\n')
        f.write('property double x\n')
        f.write('property double y\n')
        f.write('property double z\n')
        f.write('end_header\n')
        for ele in fixed_points:
            f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')

    # csv にPDSの結果出力
    with open(result_csv_path, 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(fixed_points)
    


if __name__ == '__main__':
    start = time.time()  # 時間計測用
    main() # 点群生成
    
    # 時間計測結果の表示
    elapsed_time = time.time() - start 
    print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")

