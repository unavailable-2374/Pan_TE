import time
import pickle
import logging
from pathlib import Path
from typing import Any, Callable

logger = logging.getLogger(__name__)

class RobustRunner:
    """增强的错误处理和恢复机制"""
    
    def __init__(self, config):
        self.config = config
        self.checkpoint_dir = Path(config.checkpoint_dir)
        self.checkpoint_dir.mkdir(exist_ok=True)
    
    def run_with_retry(self, func: Callable, max_retries: int = None, 
                      backoff: int = None) -> Any:
        """带重试的函数执行"""
        max_retries = max_retries or self.config.max_retries
        backoff = backoff or self.config.retry_backoff
        
        last_exception = None
        for attempt in range(max_retries):
            try:
                return func()
            except Exception as e:
                last_exception = e
                if attempt < max_retries - 1:
                    wait_time = backoff ** attempt
                    logger.warning(f"Attempt {attempt + 1} failed: {e}")
                    logger.info(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    logger.error(f"All {max_retries} attempts failed")
        
        raise last_exception
    
    def run_with_checkpoint(self, func: Callable, checkpoint_name: str) -> Any:
        """带检查点的函数执行"""
        checkpoint_file = self.checkpoint_dir / f"{checkpoint_name}.pkl"
        
        # 尝试从检查点恢复
        if checkpoint_file.exists():
            logger.info(f"Restoring from checkpoint: {checkpoint_name}")
            try:
                with open(checkpoint_file, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                logger.warning(f"Failed to restore checkpoint: {e}")
        
        # 执行函数
        result = func()
        
        # 保存检查点
        try:
            with open(checkpoint_file, 'wb') as f:
                pickle.dump(result, f)
            logger.info(f"Saved checkpoint: {checkpoint_name}")
        except Exception as e:
            logger.warning(f"Failed to save checkpoint: {e}")
        
        return result
    
    def cleanup_checkpoints(self):
        """清理所有检查点文件"""
        for checkpoint_file in self.checkpoint_dir.glob("*.pkl"):
            checkpoint_file.unlink()
        logger.info("Cleaned up all checkpoints")