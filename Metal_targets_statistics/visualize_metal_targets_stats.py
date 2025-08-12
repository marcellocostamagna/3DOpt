import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from collections import defaultdict
from pathlib import Path

# ---- Set all fonts globally to Arial, size 10 ----
matplotlib.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

def read_stats_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    metals = defaultdict(int)
    geometries = defaultdict(int)
    
    current_section = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        if line == "METALS":
            current_section = "metals"
            continue
        elif line == "GEOMETRIES":
            current_section = "geometries"
            continue
        
        if current_section == "metals" or current_section == "geometries":
            parts = line.split(": ")
            if len(parts) == 2:
                if current_section == "metals":
                    metals[parts[0]] = int(parts[1])
                else:
                    geometries[parts[0]] = int(parts[1])
    
    return metals, geometries

def plot_horizontal_histogram(data, title, svg_path, color):
    sorted_items = sorted(data.items(), key=lambda x: x[1], reverse=True)
    labels, values = zip(*sorted_items)
    labels = [label.replace("_", "\u00A0") for label in labels]  
    
    y_pos = np.arange(len(labels))
    bar_width = 0.65
    
    fig, ax = plt.subplots(figsize=(6.93, 3.3))
    bars = ax.barh(y_pos, values, color=color, edgecolor='black', linewidth=0.5, height=bar_width)
    
    for bar, value in zip(bars, values):
        ax.text(bar.get_width() + max(values)*0.01, bar.get_y() + bar.get_height()/2,
                f"{value}", va='center', fontweight='bold', fontsize=9, fontname='Arial')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{label}" for label in labels], fontsize=9, fontname='Arial')
    ax.invert_yaxis()
    ax.set_xlabel("Count", fontsize=10, fontname='Arial')
    # ax.set_title(f"{title} Distribution", fontweight='bold', fontsize=11, fontname='Arial')
    ax.tick_params(axis='both', labelsize=9)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    
    plt.tight_layout(pad=0.6)
    fig.savefig(svg_path, format='svg')
    plt.close(fig)
    print(f"✅ Plot saved to {svg_path.resolve()}")

def main():
    filename = "Metal_targets_stats.txt"
    metals, geometries = read_stats_file(filename)
    
    plot_horizontal_histogram(
        metals, 
        "Metals", 
        Path("stats_metals.svg"),
        color="#1f77b4"
    )
    plot_horizontal_histogram(
        geometries, 
        "Geometries", 
        Path("stats_geometries.svg"),
        color="#ff7f0e"
    )

if __name__ == "__main__":
    main()
