#!/usr/bin/env python3
"""
便捷脚本：迁移旧的Phase1结果到新架构
"""

import logging
import sys
from pathlib import Path
from legacy_adapter import convert_legacy_phase1

def main():
    if len(sys.argv) < 2:
        print("用法:")
        print("  python migrate_legacy_phase1.py <旧phase1结果文件> [输出文件]")
        print()
        print("示例:")
        print("  python migrate_legacy_phase1.py old_phase1_results.pkl")
        print("  python migrate_legacy_phase1.py old_phase1_results.pkl converted_phase1.pkl")
        print()
        print("支持的文件格式: .pkl, .json")
        return
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    # 设置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    try:
        print(f"🔄 正在转换旧的Phase1结果: {input_file}")
        
        converted_file = convert_legacy_phase1(input_file, output_file)
        
        print(f"✅ 转换成功!")
        print(f"📁 转换后的文件: {converted_file}")
        print()
        print("现在您可以使用以下命令运行新的流水线:")
        print(f"  python main_with_legacy.py --load-phase1 {converted_file} -g genome.fa -o output_dir")
        
    except Exception as e:
        print(f"❌ 转换失败: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()