import cv2
import numpy as np

def PickCloud(img):
    kernel = np.ones((3,3),np.uint8)
    mask_op = cv2.morphologyEx(img,
            cv2.MORPH_OPEN, kernel, iterations=3)
    mask_opcl = cv2.morphologyEx(mask_op,
            cv2.MORPH_CLOSE, kernel, iterations=10)
    mask_edge = cv2.findContours(mask_opcl, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    areas = [cv2.contourArea(res) for res in mask_edge[0]]
    where_largest = areas.index(max(areas))
    largest_edge = mask_edge[0][where_largest]
    # just edge, not containing the body
    x,y,w,h = cv2.boundingRect(largest_edge)
    boxmask = np.zeros(mask_opcl.shape)
    boxmask[y:y+h,x:x+w] = 1.
    return np.logical_and(boxmask, mask_opcl)*1.


