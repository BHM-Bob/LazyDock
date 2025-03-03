'''
Date: 2025-03-02 22:43:12
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-02 22:43:15
Description: 
'''
import numpy as np


def vectorized_sliding_average(data, window_size):
    """
    沿n_frames轴进行滑动平均（向量化实现）
    
    参数：
    data : numpy.ndarray, 形状为[n_frames, m]
    window_size : int, 滑动窗口大小
    
    返回：
    smoothed : numpy.ndarray, 形状与data相同，包含滑动平均值
    """
    n_frames, m = data.shape
    if window_size < 1:
        raise ValueError("窗口大小必须至少为1")
    if window_size > n_frames:
        raise ValueError("窗口大小不能超过数据长度")
    
    # 计算累积和
    cum_sum = np.cumsum(data, axis=0)
    
    # 构造位移后的累积和（用于计算窗口总和）
    shifted_cum_sum = np.zeros_like(cum_sum)
    shifted_cum_sum[window_size:] = cum_sum[:-window_size]
    
    # 计算每个位置的窗口总和
    window_sums = cum_sum - shifted_cum_sum
    
    # 生成动态除数数组（处理前window_size-1个位置的动态窗口）
    divisors = np.minimum(np.arange(n_frames)[:, None] + 1, window_size)
    
    # 计算滑动平均
    smoothed = window_sums / divisors
    
    return smoothed