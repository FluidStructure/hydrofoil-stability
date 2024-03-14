
# https://sectionproperties.readthedocs.io/en/stable/index.html

from sectionproperties.analysis import Section
#from sectionproperties.pre.library import circular_section
from sectionproperties.pre      import Geometry, Material
from numpy                      import loadtxt,zeros,shape
from shapely                    import Polygon

#import matplotlib.pyplot as plt
#airfoil=loadtxt("0012.txt")
#materialflg=1;

N=shape(x); airfoil=zeros((N[0],2))
airfoil[:,0]=x; airfoil[:,1]=y; airfoil=Polygon(airfoil)

steel = Material(
    name="Steel",
    elastic_modulus=200e3,  # N/mm^2 (MPa)
    poissons_ratio=1/3,  # unitless
    density=7.85e-6,  # kg/mm^3
    yield_strength=500,  # N/mm^2 (MPa)
    color="green",
)

if materialflg==1:
    geom = Geometry(geom=airfoil,material=steel)
else:
    geom = Geometry(geom=airfoil)

#geom.plot_geometry(dpi=400)

geom.create_mesh(mesh_sizes=[250])
#geom.create_mesh(mesh_sizes=[25])

sec = Section(geometry=geom)
#sec.display_mesh_info()
#sec.plot_mesh(dpi=400,materials=False)

sec.calculate_geometric_properties()
sec.calculate_warping_properties()
sec.calculate_plastic_properties()

#sec.display_results()

x_se,y_se=sec.get_sc(); area=sec.get_area(); cx, cy =sec.get_c()
if materialflg==1: # material
    # https://sectionproperties.readthedocs.io/en/stable/examples/materials/composite_analysis.html
    ixx_c, iyy_c, ixy_c = sec.get_eic(e_ref=steel)
    j = sec.get_ej(e_ref=steel)
else:  # no material
    ixx_c, iyy_c, ixy_c = sec.get_ic()
    j = sec.get_j()

print(f"xc, yc = {cx:.3f}, {cy:.3f}")
print(f"area = {area:.3f}")
print(f"Ixx, Iyy = {ixx_c:.3e}, {iyy_c:.3e}")
print(f"J = {j:.3e}")
print(f"x_sc, y_sc = {x_se:.3f}, {y_se:.3f}")

#sec.plot_centroids(dpi=400)
