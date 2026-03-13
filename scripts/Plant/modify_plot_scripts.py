import os
import glob
import re

dir_path = r'E:\CausalSDMs\scripts\Plant\plot'
files = glob.glob(os.path.join(dir_path, '*'))

replacements = {
    'EcoISEA3H': 'Plant',
    'case2_eco': 'case4_plant',
    r'\banimal\b': 'plant',
    r'\bAnimal\b': 'Plant',
    r'\bmammal\b': 'plant',
    r'\bMammal\b': 'Plant',
    r'_eco\.': '_plant.',
    r'_eco\b': '_plant',
    r'\(Eco\)': '(Plant)',
    r'\bEco\b': 'Plant'
}

for filepath in files:
    filename = os.path.basename(filepath)
    if filename == 'fig_china_cast_species_maps.R':
        continue
    if not (filename.endswith('.py') or filename.endswith('.R')):
        continue
        
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
        
    new_content = content
    for old, new in replacements.items():
        if old.isalpha() or old in ['EcoISEA3H', 'case2_eco']:
            new_content = re.sub(old, new, new_content)
        else:
            new_content = re.sub(old, new, new_content)
            
    if new_content != content:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
        print(f"Updated {filename}")
    else:
        print(f"No changes for {filename}")

# Also rename any files containing '_eco' or '_eco01' to '_plant'
for filepath in glob.glob(os.path.join(dir_path, '*')):
    filename = os.path.basename(filepath)
    if '_eco01' in filename or '_eco' in filename:
        new_filename = filename.replace('_eco01', '_plant01').replace('_eco', '_plant')
        new_filepath = os.path.join(dir_path, new_filename)
        # Avoid renaming to the same name
        if filepath != new_filepath:
            if os.path.exists(new_filepath):
                os.remove(new_filepath)
            os.rename(filepath, new_filepath)
            print(f"Renamed {filename} to {new_filename}")

print("Done!")
