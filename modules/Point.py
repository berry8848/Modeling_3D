# coding: utf-8
import numpy as np

class Point:
    def __init__(self, point):
        self.point = point
        self.x = None
        self.y = None
        self.z = None
        self.node = None
        self.coordinate = None
        self.stress = None

    def system_guid_obj_to_coordinate(self):
        '''System.GuidObjectをRhino.Objectに変換。クラス変数に座標値を保持させる'''
        p = self.point
        self.node = p[0]
        self.x = p[1]
        self.y = p[2]
        self.z = p[3]
        self.stress = np.clip(p[4], None, 0) #応力値が正の値の場合、0に変換する（常に負の値をとるようにする）
        self.coordinate = [p[1], p[2], p[3]]
   

    def change_coordinate(self):
        '''外接円を構成する点を変換するメソッド'''
        p = self.point
        self.node = p[0]
        self.x = p[1]
        self.y = p[2]
        self.z = p[3]
        self.coordinate = [p[1], p[2], p[3]]
    
    def pds_coordinate(self):
        p = self.point
        self.x = p[0]
        self.y = p[1]
        self.z = p[2]
        self.coodinate = [p[0], p[1], p[2]]
