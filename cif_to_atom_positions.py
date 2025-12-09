import math, numpy as np

def mod1(x):
    v = x - math.floor(x)
    if abs(v - 1.0) < 1e-12:  
        return 0.0
    return v

ops = [
  lambda x,y,z: (x,y,z),
  lambda x,y,z: (-x+y,-x,0.5-z),
  lambda x,y,z: (x-y,x,0.5+z),
  lambda x,y,z: (y,-x+y,-z),
  lambda x,y,z: (-y,x-y,z),
  lambda x,y,z: (x,y,0.5-z),
  lambda x,y,z: (-x,-y,0.5+z),
  lambda x,y,z: (x-y,x,-z),
  lambda x,y,z: (-x+y,-x,z),
  lambda x,y,z: (-y,x-y,0.5-z),
  lambda x,y,z: (y,-x+y,0.5+z),
  lambda x,y,z: (-x,-y,-z)
]

sites = [
  ("CaI",  (0.33333,0.66667,0.00140), 1.0, "Ca"),
  ("CaII", (0.24650,0.99330,0.25000), 1.0, "Ca"),
  ("P",    (0.39850,0.36840,0.25000), 1.0, "P"),
  ("OI",   (0.32820,0.48460,0.25000), 1.0, "O"),
  ("OII",  (0.58710,0.46490,0.25000), 1.0, "O"),
  ("OIII", (0.34340,0.25790,0.07040), 1.0, "O"),
  ("O-h",  (0.00000,0.00000,0.19600), 0.5, "O"),
  ("H",    (0.00000,0.00000,0.06080), 0.5, "H"),
]

ROUND_DECIMALS = 4

unique = {}   
order = []    

for label,(x,y,z),occ,el in sites:
    for op in ops:
        xr, yr, zr = op(x,y,z)
        xf = mod1(xr)
        yf = mod1(yr)
        zf = mod1(zr)
        key = (round(xf, ROUND_DECIMALS), round(yf, ROUND_DECIMALS), round(zf, ROUND_DECIMALS))
        if key not in unique:
            frac = (xf, yf, zf)
            unique[key] = {
                "site": label,
                "element": el,
                "occupancy": occ,
                "fx": frac[0],
                "fy": frac[1],
                "fz": frac[2],
            }
            order.append(key)
            
for i, key in enumerate(order, start=1):
    d = unique[key]
    print(f"{i}  {d['site']}  {d['element']} {d['fx']:.4f}  {d['fy']:.4f}  {d['fz']:.4f}  {d['occupancy']:.1f}")

print()
print("Total unique fractional positions (rows):", len(order))
print("Sum of occupancies (atoms in cell):", sum(d["occupancy"] for d in unique.values()))
