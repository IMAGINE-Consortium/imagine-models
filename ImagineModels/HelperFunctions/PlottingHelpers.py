import numpy as np
import matplotlib.pyplot as plt

try:
    import cmasher as cm
    has_cmasher =  True
except ImportError:
    has_cmasher = False
    

plt.ion()

def plot_slice(array, slice_dim, shp, rfp, inc, vmin, vmax, vec_dim=0, show_cbar=True, show_labels=True, label=None,
               save_fig=False, field_name=None, quiver=False, amplitude=False, cut_index =None,  plot_earth=True, cmap=None):
    """
    produces a plot showing a slice through a IMAGINE model output. The slice goes throgh the middle of array in the respective dimension slice_dim. All distances are kpc.

    :param array: vector or scalar field in grid form e.g. from on_grid()
    :param slice_dim: dimension on which to slice, e.g. 0 will give a slice through the y-z plane
    :param shp: shape of the image to plot, given as a list containing number of cells in [x, y, z]
    :param rfp: reference point = location of the smallest coordinate value in each direction, as list [refx, refy, refz]
    :param inc: increment in each dimension = size of cells used for plot, as list [incx, incy, incz]
    :param vmin: minimum value on color bar (magnetic field in muG or density in 1/cm^3)
    :param vmax: maximum value on color bar (magnetic field in muG or density in 1/cm^3)
    :param vec_dim: which dimension of vec(B) to display, can be 0, 1, 2. Necessary to plot vector field components, no effect for scalar field or if amplitude=True
    :param show_cbar: show color bar
    :param show_labels: show labels on plotted axes, default=True
    :param label: label to plot, if None, default labels are used, default=None
    :param save_fig: saves figure automatically, default=False
    :param field_name: if not None: displays name of field in plot, default=None
    :param quiver: adds field vectors displayed as arrows to plot, has no effect for scalar fields, default=False
    :param plot_earth: adds a marker on Earth's position, default=False
    :param cmap: string to select a matplotlib colormap, the default
    """
   #
    is_vector = False
    if isinstance(array, list): 
        is_vector = True 
   
    if cmap is None:
        cmap = "plasma"
        if is_vector:
            if amplitude:
                cmap = 'RdBu_r'
            else:
                cmap = "PuOr_r"
    
    fig, ax = plt.subplots()
    
    slices = [slice(0, shp[0], None), slice(0, shp[1], None), slice(0, shp[2], None)]
    if cut_index is None:
        try:
            cut_index = int(shp[slice_dim]/2)
        except IndexError:
            raise IndexError('slice_dim must be 0, 1, or 2!')
    slices[slice_dim] = slice(cut_index, cut_index +1, 1)
    slices = tuple(slices)

    sliced_at = rfp[slice_dim]+cut_index*inc[slice_dim]

    dims = [0, 1, 2]
    dims.remove(slice_dim)
    dims_label = ['x', 'y', 'z']
    
    coord1 = [i for i in range(shp[dims[0]])]
    coord2 = [i for i in range(shp[dims[1]])]
    
    if is_vector:
        if amplitude:
            comp_label = 'amplitude'
            x1 = np.array(coord1)*inc[dims[0]] + rfp[dims[0]] + inc[dims[0]] / 2.
            x2 = np.array(coord2)*inc[dims[1]] + rfp[dims[1]] + inc[dims[1]] / 2.
            vec1, vec2 = np.meshgrid(x1, x2)
            vec = np.zeros((3, np.shape(vec1)[0], np.shape(vec1)[1]))
            vec[dims[0], :, :] = vec1
            vec[dims[1], :, :] = vec2
            vec[slice_dim, :, :] = sliced_at

            B1 = np.squeeze(array[0][slices])  # need Bx and By to calculate Bphi
            B2 = np.squeeze(array[1][slices])
    
            # calculate cross product between vec(position) and vec(B_horizontal) to get direction of Bphi 
            cross = np.cross(np.array([B1.T, B2.T, np.zeros_like(B1.T)]), vec, axis=0)
            sign = np.sign(cross[-1, :, :])
            array_plot = np.squeeze(np.linalg.norm(array, axis=0)[slices]).T*np.squeeze(sign)
        else:
            comp_label = dims_label[vec_dim] + '_component'
            array_plot = np.squeeze(array[vec_dim][slices]).T
    else:
        array_plot = np.squeeze(array[slices]).T
        
    ax.imshow(array_plot, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)

    if is_vector:    
        if quiver:  # add arrows indicating B-field direction
            # only every fifth cell gets an arrow for better overview
            B1_quiv = np.squeeze(array[dims[0]][slices])
            B2_quiv = np.squeeze(array[dims[1]][slices])
            plt.quiver(coord1[::5], coord2[::5], B1_quiv[::5, ::5].T, B2_quiv[::5, ::5].T, pivot='mid')
    
    slice_dim_label = dims_label[slice_dim]
    dims_label.remove(slice_dim_label)

    if plot_earth:
        earth = [-8.5, 0, 0]
        plt.scatter((abs(rfp[dims[0]]-earth[dims[0]])) / (inc[dims[0]]*shp[dims[0]]) * shp[dims[0]],
                    (abs(rfp[dims[1]]-earth[dims[1]])) / (inc[dims[1]]*shp[dims[1]]) * shp[dims[1]],
                    marker='o', s=20, c='0.5')
    
    xticks_label = ['%.2f' % (rfp[dims[0]] + i*shp[dims[0]]/4*inc[dims[0]]) for i in range(5)]
    if float(int(rfp[dims[0]])) != 0:
        if rfp[dims[0]] / float(int(rfp[dims[0]])) - 1 < 1e-3:
            xticks_label = [int(rfp[dims[0]] + i*shp[dims[0]]/4*inc[dims[0]]) for i in range(5)]
    xticks_loc = [i*shp[dims[0]]/4 for i in range(5)]

    yticks_label = ['%.2f' % (rfp[dims[1]] + i*shp[dims[1]]/4*inc[dims[1]]) for i in range(5)]
    if float(int(rfp[dims[1]])) != 0:
        if rfp[dims[1]] / float(int(rfp[dims[1]])) - 1 < 1e-3:
            yticks_label = [int(rfp[dims[1]] + i*shp[dims[1]]/4*inc[dims[1]]) for i in range(5)]
    yticks_loc = [i*shp[dims[1]]/4 for i in range(5)]
        

    if show_labels:
        ax.set_yticks(yticks_loc, labels=yticks_label)
        ax.set_ylabel(r"$%s$ / kpc" % dims_label[1])
        ax.set_xticks(xticks_loc, labels=xticks_label)
        ax.set_xlabel(r"$%s$ / kpc" % dims_label[0])
        keyword = '%s = %.3f kpc' % (slice_dim_label, sliced_at) 
        keyword = field_name + ', ' + keyword if field_name is not None else keyword
        ax.text(0, 1, keyword, va='bottom', transform = ax.transAxes)
    else:
        plt.axis('off')


    if show_cbar:
        cbar = fig.colorbar(ax.images[0], orientation="horizontal", shrink=0.5, aspect=30)
        if show_labels:
            if label is None:
                if is_vector:
                    if amplitude:
                        label = r'$|B|$ sign($B_\phi$) / $\mu$G'
                    else:
                        label = r'$B_%s$ / $\mu$G' % comp_label[0]
                else: 
                    label = r'$n_\mathrm{th}$ / $\mathrm{cm}^-3$'                 
            cbar.set_label(label)
    plt.tight_layout()
    if save_fig:
        name = comp_label + "_" + dims_label[0] + "_" + dims_label[1] + "_plane"
        if not is_vector:
            name = name[1:]
        plt.savefig(name)    
    plt.show()
    plt.close()