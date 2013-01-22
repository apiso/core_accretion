"""
read in data generated by shoot2.py for plotting purposes.
"""

id = "mc5a10_500"

getdat = 1
if getdat:
    from utils.constants import Me, Re
    from RadSGPoly import shoot2 as s2
    datfolder = s2.datfolder
    file = id + ".npz"
    params, arr = s2.parammodelsload(file)
    Matm = arr.M - params.mco
