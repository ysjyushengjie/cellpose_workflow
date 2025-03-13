"""
======================================================
@File   : 02_based_cellpose_segment_cell.py
@Author : ouyangkang
@Date   : 2025.03.04
@Email  : ouyangkang@genomics.cn
@Desc   : utilizing cellpose model to segment cell
@Version: v1
@Update : 2025.03.13
======================================================
update record：
- 2025.03.13：adding annotation information and modifying error
"""

from cellpose import models, io
import matplotlib.pyplot as plt
import sys


def cell_segment(registered_file, chip_id, model_path='/data/work/ascites/model/models/cyto3'):
    """
    utilizing cellpose model cyto3 to segment cell in HE image
    
    Parameters:
        registered_file: the absolute pathway of registered HE image(str)
        chip_id: the chip id (str)
        model_path: the path of cellpose model
        
    Return:
        without return
    """
    
    # read registered file
    img = io.imread(registered_file)
    # load pre-trained model
    model = models.CellposeModel(gpu=True, pretrained_model=model_path)
    # run pre-trained model
    masks, flows, styles = model.eval(img, diameter=None, channels=[0,3])
    
    io.masks_flows_to_seg(img, masks, flows, registered_file, styles, channels=[0,3])
    io.save_masks(img, masks, flows, f'{chip_id}_cellpose.tif', png=True, tif=True, channels=[0,3])


if __name__ == "__main__":
    # regesited tiff, for example '/data/work/ascites/aligned_picture/Y00722A9_HE_aligned_OYK.tif'
    input_file = sys.argv[1]
    model_path = sys.argv[2]
    chip_id = sys.argv[3]
    
    cell_segment(input_file, chip_id, model_path)