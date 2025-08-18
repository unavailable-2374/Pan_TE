import hashlib
import pickle
import logging
from pathlib import Path
from functools import wraps

logger = logging.getLogger(__name__)

def cache_result(func):
    """缓存计算密集型操作结果的装饰器"""
    cache_dir = Path("cache")
    cache_dir.mkdir(exist_ok=True)
    
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        # 生成缓存键
        cache_key = hashlib.md5(
            f"{func.__name__}_{str(args)}_{str(kwargs)}".encode()
        ).hexdigest()
        cache_file = cache_dir / f"{func.__name__}_{cache_key}.pkl"
        
        # 检查缓存
        if cache_file.exists():
            try:
                logger.debug(f"Loading cached result for {func.__name__}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                logger.warning(f"Failed to load cache: {e}")
        
        # 执行函数
        result = func(self, *args, **kwargs)
        
        # 保存缓存
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(result, f)
            logger.debug(f"Cached result for {func.__name__}")
        except Exception as e:
            logger.warning(f"Failed to save cache: {e}")
        
        return result
    
    return wrapper

def clear_cache():
    """清理所有缓存文件"""
    cache_dir = Path("cache")
    if cache_dir.exists():
        for cache_file in cache_dir.glob("*.pkl"):
            cache_file.unlink()
        logger.info("Cleared all cache files")