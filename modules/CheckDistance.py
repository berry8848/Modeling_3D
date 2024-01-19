import math

class CheckDistance:
    '''応力値から点間距離を求め，点間距離内に他の点が含まれるか否かを判断する'''
    def __init__(self, ALLOWABLE_STRESS, RADIUS):
        self.ALLOWABLE_STRESS = ALLOWABLE_STRESS
        self.r = RADIUS

    def check_distance(self, fixed_points, candidate_point, stress):
        density = stress_to_density(stress, self.ALLOWABLE_STRESS) # 応力から密度の変換
        long = density_to_long(density, self.r) # 密度から点間距離の変換
        # 点間距離内に他の点が含まれているか否かを判定．
        # 点間距離内に他の点が含まれていたらFalseを返す
        b_x = candidate_point[0]
        b_y = candidate_point[1]
        b_z = candidate_point[2]
        check = True
        for a in fixed_points:
            a_x = a[0]
            a_y = a[1]
            a_z = a[2]
            distance = (a_x- b_x)**2 + (a_y- b_y)**2 + (a_z- b_z)**2
            if long**2 > distance:
                check = False
                break
        return check
        
# 応力と密度の関係式．※関係式が微妙なため，臨時で別の関数
def stress_to_density(stress, ALLOWABLE_STRESS):
    density = abs(stress) / ALLOWABLE_STRESS
    return density

# 密度と点間距離の関係式
def density_to_long(density, r):
    long = 6 * r * math.sqrt(math.sqrt(2)*math.pi/ density)
    return long
