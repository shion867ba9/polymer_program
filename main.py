import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path = 'C:\Users\shion\PycharmProjects\polymerization_program'

# === 初期パラメーター ===

meoh_ratio = 0.25    # メタ比
temp = 60    # 重合温度
target_conversion = 30    # 目標重合度
div = 100
conversion = [x/div/100 for x in range(0, target_conversion*(div+1), target_conversion)]
dist = target_conversion/div/100
# --- 仮パラメータ ---
x = 0
T = 0
Pv = 0


# === 開始剤 ===
initiator = 'AIBN'
f = 0.6    # 開始剤効率

# === 連鎖移動定数 ===
Ctr = 7.7/np.power(10, 4)
Cm = 1.35/np.power(10, 4)
Cs = 4.3/np.power(10, 4)
Ci = 0

# === 速度定数 ===
kd = 0.85/np.power(10, 5)
kp = 3700/np.power(10, 4)
kt = 58/np.power(10, 6)
ktr = 7.7/np.power(10, 4)*kp


# === 分子量 ===
fw_meoh = 32.04
fw_vac = 86.09
fw_pvac = 86.09
fw_init = 121    # 開始剤(AIBN)の分子量

# === density ===
d_meoh = 0.7
d_vac = 0.93
d_pvac = 1.2

Pv_list = []
T_list = []

for x in conversion:
    # --- unit: g ---
    amount_whole = 1600    # 全量
    amount_vac_0 = amount_whole*(1-meoh_ratio)    # 酢ビ重量
    amount_vac = lambda x: amount_vac_0*(1-x)
    amount_pvac = lambda x: amount_vac_0*x
    amount_meoh_0 = amount_whole*meoh_ratio    # メタノール重量
    # amount_meoh_feed
    amount_meoh = amount_meoh_0 # + amount_meoh_feed
    amount_init_0 = 1


    # --- unit: mol ---
    as_meoh = amount_meoh/fw_meoh
    as_vac = amount_vac(x)/fw_vac
    as_pvac = amount_pvac(x)/fw_pvac
    as_init = amount_init_0/fw_init

    # --- unit: L ---
    vol_meoh = amount_meoh/d_meoh
    vol_vac = amount_vac(x)/d_vac
    vol_pvac = amount_pvac(x)/d_pvac
    vol_total = (vol_meoh + vol_vac + vol_pvac)/1000

    # --- concentration: mol/L ---
    conc_meoh = as_meoh / vol_total
    conc_vac = as_vac / vol_total
    conc_pvac = as_pvac / vol_total
    conc_init = as_init / vol_total

    Rp = kp*conc_vac*np.power(kd*conc_init*f/kt, 1/2)
    dT = as_vac*dist/Rp
    T += dT
    delta_inv_Pn = Cm + Ci*conc_init/conc_vac + Cs*conc_meoh/conc_vac # + Ctr*conc_trans/conc_vac
    delta_Pv = 1/delta_inv_Pn * 1.824
    a = x / (x + dist)
    Pv = np.power(a*np.power(Pv, 0.62)+(1-a)*np.power(delta_Pv, 0.62), 1/0.62)
    inv_Pv = 1/Pv
    T_list.append(T)
    Pv_list.append(Pv)

print(delta_Pv, Pv, T)
plt.plot(T_list, Pv_list)
