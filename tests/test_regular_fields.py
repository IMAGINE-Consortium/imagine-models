import ImagineModels as img

import pytest
import numpy as np


regular_models = ['JaffeMagneticField', 'HelixMagneticField', 'JF12RegularField', 'SunMagneticField'
                  #'AxiSymmetricSpiral'
                  ]

known_positions = {'JF12RegularField': [([0, 0, 0], [0, 0, 0]), # Galactic center
                                        ([.1, .3, .4], [0, 0, 0]), # Within inner boundary
                                        ([3, 20., .3], [0, 0, 0, ]),], # outside outer boundary 
                   'JaffeMagneticField': [([0, 0, 0], [0., 0., 0.]),], # origin 
                   # 'AxiSymmetricSpiral': [([0, 0, 0], [0, 0, 0]), # origin
                   #                       ], 
                   'HelixMagneticField': [([0, 0, 0], [0, 0, 0]), # origin
                                         ],
                    'SunMagneticField': [([0, 0, 0], [0.41582338163551863, 1.9562952014676114, 0.0]), # origin
                                         ],
                   }

# irregular_grid
grid_x = np.asarray([-2.2, -1., 0., .1])
grid_y = np.asarray([-1., 0., 3.])
grid_z = np.asarray([-5., 0., .1])

shape = [2, 3, 4]
increment = [.1, 3., .1]
zeropoint = [-2., 3., .1]

def test_at_position():
    for model_string in regular_models:
        mo = getattr(img, model_string) 
        mo_1 = mo()
        mo_2 = mo(grid_x, grid_y, grid_z)
        mo_3 = mo(shape,  zeropoint, increment)
        
        for position, value in known_positions[model_string]:
            x_pos, y_pos, z_pos = position
            pos_1 = mo_1.at_position(x_pos, y_pos, z_pos)
            pos_2 = mo_2.at_position(x_pos, y_pos, z_pos)
            pos_3 = mo_3.at_position(x_pos, y_pos, z_pos)
            for j in range(3):
                assert pos_1[j] == value[j]
                assert pos_2[j] == value[j]
                assert pos_3[j] == value[j]
            
            
def test_grid():
    # this function tests the consistency of the different ways to initialize a model
    for model_string in regular_models:
        mo = getattr(img, model_string) 
        mo_1 = mo()
        mo_2 = mo(grid_x, grid_y, grid_z)
        mo_3 = mo(shape, zeropoint, increment)
        
        with pytest.raises(ValueError):
            mo_1.on_grid()
        b_field_1 = mo_1.on_grid(grid_x, grid_y, grid_z)
        b_field_2a = mo_2.on_grid()
        b_field_2b = mo_2.on_grid(shape=shape, reference_point=zeropoint, increment=increment)
        b_field_3 = mo_3.on_grid()
        for j in range(3):
            assert np.array_equal(b_field_1[j], b_field_2a[j], equal_nan=True)
            assert np.array_equal(b_field_2b[j], b_field_3[j], equal_nan=True)
            
        
def test_interface():
    grid_x = np.asarray([-2.2, -1., 0., .1])
    grid_y = np.asarray([-1., 0., 3.])
    grid_z = np.asarray([-5., 0., .1])

    shape = [2, 3, 4]
    increment = [.1, 3., .1]
    zeropoint = [-2., 3., .1]

    umf = img.UniformMagneticField()

    umf.on_grid(grid_x, grid_y, grid_z)
    umf.on_grid(shape=shape, reference_point=zeropoint, increment=increment)
    umf.on_grid(shape=np.asarray(shape), reference_point=zeropoint, increment=increment)
    
    float_shape = [2., 3., 4.]
    one_int_inc = [2.1, 3, .1]
    one_float_shape = [2, 3, 4.]
    
    with pytest.raises(TypeError):
        umf.on_grid(shape, zeropoint, increment) # must be passed by keyword
    with pytest.raises(TypeError):
        umf.on_grid(shape=float_shape, reference_point=zeropoint, increment=increment) # shape must be int
    with pytest.raises(TypeError):
        umf.on_grid(shape=one_float_shape, reference_point=zeropoint, increment=increment) # all shape must be int
    with pytest.raises(TypeError):
        umf.on_grid(shape=shape, reference_point=zeropoint, increment=one_int_inc) # all increments must be float
    with pytest.raises(TypeError):
        umf.on_grid(grid_x.tolist(), grid_y, grid_z) # lists won't work
    

            

def test_parameter_update():
    umf = img.UniformMagneticField()
    assert umf.bx == 0.
    assert umf.by == 0.
    assert umf.bz == 0.
    
    assert umf.at_position(2.4, 2.1, -.2) == (0., 0., 0.)
    
    umf.bx = -3.2 
    
    assert umf.bx == -3.2 
    
    assert umf.at_position(2.4, 2.1, -.2) == (-3.2, 0., 0.)