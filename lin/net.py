#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

# Example transition matrix
data = {
    "GeneA": [0.0, 0.5, 0.2, 0.1, 0.05, 0.0, 0.0],
    "GeneB": [0.5, 0.0, 0.5, 0.2, 0.1, 0.1, -0.3],
    "GeneC": [0.2, 0.5, 0.0, 0.4, -0.2, -0.1, 0.3],
    "GeneD": [0.1, 0.2, 0.4, 0.0, 0.2, 0.2, 0.0],
    "GeneE": [0.05, 0.1, -0.2, 0.2, 0.0, 0.3, 0.25],
    "GeneF": [0.0, 0.1, -0.1, 0.2, 0.3, 0.0, 0.0],
    "GeneG": [0.0, -0.3, 0.3, 0.0, 0.25, 0.0, 0.0]
}

transition_df = pd.DataFrame(data, index=["GeneA", "GeneB", "GeneC", "GeneD","GeneE","GeneF","GeneG"])

G = nx.DiGraph()

# Add edges with weights (probabilities)
for from_gene in transition_df.index:
    for to_gene in transition_df.columns:
        weight = transition_df.loc[from_gene, to_gene]
        if abs(weight) > 0:
            G.add_edge(from_gene, to_gene, weight=weight)

# Define layout (spring_layout or circular_layout)
pos = nx.spring_layout(G, seed=42)

# Get edge weights
weights = [G[u][v]['weight'] * 5 for u, v in G.edges()]  # scale up for visibility

# #### option 1: thickness represents weight
# plt.figure(figsize=(8, 6))
# nx.draw(
#     G,
#     pos,
#     with_labels=True,
#     node_size=2000,
#     node_color="lightblue",
#     edge_color="gray",
#     width=weights,
#     font_size=12,
#     arrows=True,
# )
# edge_labels = {(u, v): f"{G[u][v]['weight']:.2f}" for u, v in G.edges()}
# nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="red")

# plt.title("Gene Interaction Network (Edge Thickness = Transition Probability)")
# plt.show()

#### option 2: color represents weight
# Choose a color map ('plasma', 'coolwarm', 'viridis', 'inferno', etc.)
cmap = plt.cm.coolwarm

# Normalize weights between 0 and 1 for colormap scaling
norm = plt.Normalize(vmin=min(weights), vmax=max(weights))
edge_colors = [cmap(norm(w)) for w in weights]

fig, ax = plt.subplots(figsize=(12, 10), constrained_layout=True)

nx.draw(
    G,
    pos,
    with_labels=True,
    node_size=2000,
    node_color="lightgrey",
    edge_color=edge_colors,
    width=2,
    font_size=12,
    arrows=True,
    ax=ax
)

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.04)
cbar.set_label("Transition Probability / Correlation Strength")

edge_labels = {(u, v): f"{G[u][v]['weight']:.2f}" for u, v in G.edges()}
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="black", font_size=10, ax=ax)

plt.title("Gene Interaction Network (Edge Color = Transition Strength)")
plt.show()

