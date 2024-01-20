# Delaunay分割後、メッシュ生成
import numpy as np
from scipy.spatial import Delaunay
import time
from modules import CrossingNumberAlgorithm
import trimesh

def create_stl_with_edges(vertices, faces, edges, stl_path):
    # vertices: 頂点座標の配列 (N x 3)
    # faces: 頂点インデックスの配列 (M x 3)
    # edges: エッジの配列 (K x 2)
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, edges=edges)
    mesh.export(stl_path, file_type='stl', include_edges=True)

start = time.time()  # 時間計測用
edges = [] # PLYファイルのedge用


Input_file = 'Output/PDS/pds3/pds3.csv' # Inputファイル
mesh_data = 'Input/Mesh_Data/Hourgrass2_mesh.stl' # 物体の表面形状データ。
Output_file = 'Output/Delaunay/delaunay3.ply' # Outputファイル

# ファイルの読み込み。物体の頂点を定義する
vertices = np.loadtxt(Input_file, delimiter=',')
print('vertices', vertices[:3])


# Delaunay分割の作成
tri = Delaunay(vertices)
print('simplices：', tri.simplices)

# CrossNumberAlgorithm
CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm(mesh_data)

print('len(tri.simplice)：', len(tri.simplices))
i = 0

# edgeノードの作成
for simplice in tri.simplices:
    flg01 = CNA.majority_vote((vertices[simplice[0]]+vertices[simplice[1]])/2)
    flg02 = CNA.majority_vote((vertices[simplice[0]]+vertices[simplice[2]])/2)
    flg03 = CNA.majority_vote((vertices[simplice[0]]+vertices[simplice[3]])/2)
    flg12 = CNA.majority_vote((vertices[simplice[1]]+vertices[simplice[2]])/2)
    flg13 = CNA.majority_vote((vertices[simplice[1]]+vertices[simplice[3]])/2)
    flg23 = CNA.majority_vote((vertices[simplice[2]]+vertices[simplice[3]])/2)

    # Trueのとき，edge生成
    if flg01: edges.append([simplice[0], simplice[1]])
    if flg02: edges.append([simplice[0], simplice[2]])
    if flg03: edges.append([simplice[0], simplice[3]])
    if flg12: edges.append([simplice[1], simplice[2]])
    if flg13: edges.append([simplice[1], simplice[3]])
    if flg23: edges.append([simplice[2], simplice[3]])

    print(i)
    i+=1



print('edge 重複削除前個数：', len(edges))
#重複した座標を削除
edges = np.unique(edges, axis=0)
print('edge 重複削除後個数：', len(edges))

# plyで保存
with open(Output_file, 'w', newline="") as f:
    f.write('ply\n')
    f.write('format ascii 1.0\n')
    f.write('comment VCGLIB generated\n')
    f.write('element vertex '+str(len(vertices))+'\n')
    f.write('property float x\n')
    f.write('property float y\n')
    f.write('property float z\n')
    f.write('element edge '+str(len(edges))+'\n')
    f.write('property int vertex1\n')
    f.write('property int vertex2\n')
    f.write('end_header\n')

    for ele in vertices: # 点群の座標値入力
        f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')
    
    for ele in edges: # 三角形を構成する点群のノード入力
        f.write(str(ele[0])+' '+str(ele[1])+'\n')

# create_stl_with_edges(vertices, None, edges, 'Output/Delaunay_Python/delaunay_mesh.stl')

# 計測結果
elapsed_time = time.time() - start
print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")

