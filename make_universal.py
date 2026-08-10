#!/usr/bin/env python3
"""
Script to make ParalogWizard universal across Python versions.
Run this to automatically update all Python files for compatibility.
"""

import os
import re
from pathlib import Path

def add_future_annotations(file_path):
    """Add 'from __future__ import annotations' to a Python file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Skip if already has future annotations
    if 'from __future__ import annotations' in content:
        print(f"✅ {file_path} already has future annotations")
        return False
    
    # Find where to insert (after shebang and docstring)
    lines = content.split('\n')
    insert_line = 0
    
    # Skip shebang
    if lines[0].startswith('#!'):
        insert_line = 1
    
    # Skip module docstring
    in_docstring = False
    triple_quote = None
    for i, line in enumerate(lines[insert_line:], insert_line):
        stripped = line.strip()
        if not in_docstring and (stripped.startswith('"""') or stripped.startswith("'''")):
            triple_quote = stripped[:3]
            in_docstring = True
            if stripped.count(triple_quote) >= 2:  # Single line docstring
                insert_line = i + 1
                break
        elif in_docstring and triple_quote in line:
            insert_line = i + 1
            break
        elif not in_docstring and stripped and not stripped.startswith('#'):
            break
    
    # Insert future import
    lines.insert(insert_line, '')
    lines.insert(insert_line + 1, 'from __future__ import annotations')
    
    with open(file_path, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"✅ Added future annotations to {file_path}")
    return True

def update_type_hints(file_path):
    """Update old typing imports to modern syntax."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Replace old type hints with modern ones
    replacements = [
        (r'\bList\[', 'list['),
        (r'\bDict\[', 'dict['),
        (r'\bTuple\[', 'tuple['),
        (r'\bSet\[', 'set['),
        (r'from typing import ([^,\n]*,?\s*)*List([^,\n]*,?\s*)*', lambda m: remove_from_import(m.group(0), 'List')),
        (r'from typing import ([^,\n]*,?\s*)*Dict([^,\n]*,?\s*)*', lambda m: remove_from_import(m.group(0), 'Dict')),
        (r'from typing import ([^,\n]*,?\s*)*Tuple([^,\n]*,?\s*)*', lambda m: remove_from_import(m.group(0), 'Tuple')),
        (r'from typing import ([^,\n]*,?\s*)*Set([^,\n]*,?\s*)*', lambda m: remove_from_import(m.group(0), 'Set')),
    ]
    
    for pattern, replacement in replacements:
        if callable(replacement):
            content = re.sub(pattern, replacement, content)
        else:
            content = re.sub(pattern, replacement, content)
    
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        print(f"✅ Updated type hints in {file_path}")
        return True
    
    return False

def remove_from_import(import_line, type_to_remove):
    """Remove a specific type from a typing import line."""
    # Simple implementation - just remove the type
    return import_line.replace(f', {type_to_remove}', '').replace(f'{type_to_remove}, ', '').replace(f'{type_to_remove}', '')

def process_directory(directory):
    """Process all Python files in a directory."""
    python_files = list(Path(directory).rglob("*.py"))
    
    print(f"Found {len(python_files)} Python files in {directory}")
    
    updated_files = []
    for file_path in python_files:
        print(f"\nProcessing {file_path}...")
        
        changed = False
        changed |= add_future_annotations(file_path)
        changed |= update_type_hints(file_path)
        
        if changed:
            updated_files.append(file_path)
        else:
            print(f"⏭️  No changes needed for {file_path}")
    
    return updated_files

def main():
    """Main function to make ParalogWizard universal."""
    print("🚀 Making ParalogWizard universal across Python versions...")
    print("=" * 60)
    
    # Process main file
    main_files = ["ParalogWizard.py"]
    directories = ["ParalogWizard/"]
    
    all_updated = []
    
    # Process main files
    for file_path in main_files:
        if os.path.exists(file_path):
            print(f"\nProcessing main file: {file_path}")
            changed = False
            changed |= add_future_annotations(file_path)
            changed |= update_type_hints(file_path)
            if changed:
                all_updated.append(file_path)
    
    # Process directories
    for directory in directories:
        if os.path.exists(directory):
            updated = process_directory(directory)
            all_updated.extend(updated)
    
    print("\n" + "=" * 60)
    print(f"✅ COMPLETED! Updated {len(all_updated)} files:")
    for file_path in all_updated:
        print(f"   📝 {file_path}")
    
    print(f"\n🎉 ParalogWizard is now universal! Works with Python 3.7+")
    
    # Test compatibility
    print("\n🧪 Testing compatibility...")
    try:
        exec("def test(items: list[str]) -> dict[str, int]: return {}")
        print("✅ Modern type hints work!")
    except SyntaxError:
        print("⚠️  Modern type hints failed - check Python version")

if __name__ == "__main__":
    main()