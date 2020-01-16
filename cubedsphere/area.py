import numpy as np
from .grid import csgrid_GMAO

def calc_cs_face_area(lon_b, lat_b, r_sphere = 6.375e6):
    """Calculate area of cubed-sphere grid cells on one face
    Inputs must be in degrees. Edge arrays must be
    shaped [N+1 x N+1]
    """
    
    # Convert inputs to radians
    lon_b_rad = lon_b * np.pi / 180.0
    lat_b_rad = lat_b * np.pi / 180.0
    
    r_sq = r_sphere * r_sphere
    n_cs = lon_b.shape[1] - 1
    
    # Allocate output array
    cs_area = np.zeros((n_cs,n_cs))
    
    # Ordering
    valid_combo = np.array([[1,2,4],[2,3,1],[3,2,4],[4,1,3]]) - 1
    
    for i_lon in range(n_cs):
        for i_lat in range(n_cs):
            lon_corner = np.zeros(4)
            lat_corner = np.zeros(4)
            xyz_corner = np.zeros((4,3))
            for i_vert in range(4):
                x_lon = i_lon + (i_vert > 1)
                x_lat = i_lat + (i_vert == 0 or i_vert == 3)
                lon_corner[i_vert] = lon_b_rad[x_lon,x_lat]
                lat_corner[i_vert] = lat_b_rad[x_lon,x_lat]
            for i_vert in range(4):
                xyz_corner[i_vert,:] = ll2xyz(lon_corner[i_vert],lat_corner[i_vert])
            tot_ang = 0.0
            for i_corner in range(4):
                curr_combo = valid_combo[i_corner,:]
                xyz_mini = np.zeros((3,3))
                for i_mini in range(3):
                    xyz_mini[i_mini,:] = xyz_corner[curr_combo[i_mini],:]
                curr_ang = sphere_angle(xyz_mini[0,:],xyz_mini[1,:],xyz_mini[2,:])
                tot_ang += curr_ang
            cs_area[i_lon,i_lat] = r_sq * (tot_ang - (2.0*np.pi))
    
    return cs_area

def ll2xyz(lon_pt,lat_pt):
    """Converts a lon/lat pair (in radians) to cartesian co-ordinates
    Vector should point to the surface of the unit sphere"""

    xPt = np.cos(lat_pt) * np.cos(lon_pt)
    yPt = np.cos(lat_pt) * np.sin(lon_pt)
    zPt = np.sin(lat_pt)
    return [xPt,yPt,zPt]

def sphere_angle(e1,e2,e3):
    # e1: Mid-point
    # e2 and e3 to either side
    pVec = np.ones(3)
    qVec = np.ones(3)
    pVec[0] = e1[1]*e2[2] - e1[2]*e2[1]
    pVec[1] = e1[2]*e2[0] - e1[0]*e2[2]
    pVec[2] = e1[0]*e2[1] - e1[1]*e2[0]

    qVec[0] = e1[1]*e3[2] - e1[2]*e3[1]
    qVec[1] = e1[2]*e3[0] - e1[0]*e3[2]
    qVec[2] = e1[0]*e3[1] - e1[1]*e3[0]
    ddd = np.sum(pVec*pVec) * np.sum(qVec*qVec)
    if ddd <= 0.0:
        angle = 0.0;
    else:
        ddd = np.sum(pVec*qVec)/np.sqrt(ddd);
        if (np.abs(ddd)>1.0):
            angle = np.pi/2.0;
        else:
            angle = np.arccos(ddd);

    return angle

def calc_cs_area(cs_grid=None,cs_res=None):
    """Return area in m2 for each cell in a cubed-sphere grid
    Uses GMAO indexing convention (6xNxN)
    """
    # Calculate area on a cubed sphere
    if cs_res is None:
        cs_res = cs_grid['lon_b'].shape[-1] - 1
    elif cs_grid is None:
        cs_grid = csgrid_GMAO(cs_res)
    elif cs_grid is not None and cs_res is not None:
        assert cs_res == cs_grid['lon_b'].shape[-1], 'Routine calc_cs_area received inconsistent inputs' 
    cs_area = np.zeros((6,cs_res,cs_res))
    cs_area[0,:,:] = calc_cs_face_area(cs_grid['lon_b'][0,:,:],cs_grid['lat_b'][0,:,:])
    for i_face in range(1,6):
        cs_area[i_face,:,:] = cs_area[0,:,:].copy()
    return cs_area
