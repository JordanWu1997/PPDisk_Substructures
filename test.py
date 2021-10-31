from ALMA_data_manipulation import ObsData
from azimuthal_PV_diagram import generate_azimuthal_PV_diagram
from generate_Keplerian_disk_PV import *

# Load ALMA data
IM_Lup = ObsData(
    '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO.fits', 158.)
IM_Lup.stellar_property(1.12, 4250)
IM_Lup.disk_property(47.5,
                     144.5,
                     np.array([117.]),
                     np.array([117. * 0.13]),
                     offset_x=-1.5,
                     offset_y=1.0)
IM_Lup.planet_property(np.array([69.]), np.array([0.77]))

# AzimuthalPVDiagram
r_center = IM_Lup.asec2pix(IM_Lup.planet_sep[0])
r_width = 1
aspect_ratio = 0.5

plt.figure(figsize=(18, 8))

generate_azimuthal_PV_diagram(IM_Lup.data,
                              IM_Lup.freqax2velax(),
                              r_center,
                              r_width,
                              IM_Lup.disk_pos,
                              IM_Lup.disk_inc,
                              IM_Lup.offset_x_pix,
                              IM_Lup.offset_y_pix,
                              new_figure=False,
                              show_plot=False)

plt.axvline(IM_Lup.planet_pos[0], label='planet PA')

plot_LOS_Kepleran_PV_cut(IM_Lup,
                         r_center,
                         361,
                         aspect_ratio,
                         new_figure=False,
                         show_plot=False)



plt.title(
    'Azimuthal Position-Velocity Diagram\n(r_center={:.1f} au ({:d} pix), r_width={:.1f} au ({:d} pix), z/r={:.2f}, star={:.1f} solar_mass)'
    .format(IM_Lup.pix2au(r_center), r_center, IM_Lup.pix2au(r_width), r_width,
            aspect_ratio, IM_Lup.stellar_mass))
plt.show()
