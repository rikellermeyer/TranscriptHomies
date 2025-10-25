#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

####################
#  Random Walk  
####################
def gene_random_walk(transition_matrix, start_gene, steps=1000):
    """
    Simulates a biased random walk between genes using a transition probability matrix.

    Args:
        transition_matrix (pd.DataFrame): Square matrix with probabilities of moving between genes.
        start_gene (str): The starting gene.
        steps (int): Number of steps to simulate.

    Returns:
        list: The sequence of visited genes.
    """
    genes = list(transition_matrix.columns)
    if start_gene not in genes:
        raise ValueError(f"Start gene {start_gene} not found in transition matrix.")

    current_gene = start_gene
    path = [current_gene]

    for _ in range(steps - 1):
        # Get probabilities for the next move
        probs = transition_matrix.loc[current_gene].values
        # Normalize in case of rounding issues
        probs = probs / probs.sum()
        # Choose next gene based on bias probabilities
        next_gene = np.random.choice(genes, p=probs)
        path.append(next_gene)
        current_gene = next_gene

    return path


# Example transition matrix
data = {
    "GeneA": [0.0, 0.5, 0.2, 0.1],
    "GeneB": [0.5, 0.0, 0.5, 0.2],
    "GeneC": [0.2, 0.5, 0.0, 0.4],
    "GeneD": [0.1, 0.2, 0.4, 0.0]
}
transition_df = pd.DataFrame(data, index=["GeneA", "GeneB", "GeneC", "GeneD"])

# Walk
path = gene_random_walk(transition_df, start_gene="GeneA", steps=1000)

# Count visits to each gene
visit_counts = pd.Series(path).value_counts(normalize=True)

# Plot the visit frequencies
plt.figure(figsize=(6,4))
visit_counts.plot(kind="bar", color="skyblue")
plt.title("Gene Visit Frequencies During Random Walk")
plt.xlabel("Gene")
plt.ylabel("Visit Frequency")
plt.grid(axis='y')
plt.show()

# Count how many times each transition occur
edges = list(zip(path[:-1], path[1:]))
edge_counts = pd.Series(edges).value_counts()


G = nx.DiGraph()
# edges with weights
for (src, dst), weight in edge_counts.items():
    G.add_edge(src, dst, weight=weight)

plt.figure(figsize=(8, 6))
pos = nx.spring_layout(G, seed=42, k=0.8)  # Layout for better spacing

# Get edge weights for thickness
weights = [G[u][v]['weight'] for u, v in G.edges()]
max_w = max(weights)
scaled_weights = [3 * (w / max_w) for w in weights]  # Scale for visibility

# nodes and edges
nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=1500)
nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')

nx.draw_networkx_edges(
    G, pos,
    width=scaled_weights,
    edge_color='gray',
    arrows=True,
    arrowsize=20,
    alpha=0.7
)

# Edge labels 
edge_labels = {edge: f"{G[edge[0]][edge[1]]['weight']}" for edge in G.edges()}
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

plt.title("Gene-Gene Transition Network from Biased Random Walk", fontsize=14)
plt.axis('off')
plt.tight_layout()
plt.show()




#################
#  Network graph 
#################
# Example transition matrix
data = {
    "GeneA": [0.0, 0.5, 0.2, 0.1],
    "GeneB": [0.5, 0.0, 0.5, 0.2],
    "GeneC": [0.2, 0.5, 0.0, 0.4],
    "GeneD": [0.1, 0.2, 0.4, 0.0]
}
transition_df = pd.DataFrame(data, index=["GeneA", "GeneB", "GeneC", "GeneD"])

G = nx.DiGraph()

# Add edges with weights (probabilities)
for from_gene in transition_df.index:
    for to_gene in transition_df.columns:
        weight = transition_df.loc[from_gene, to_gene]
        if weight > 0:
            G.add_edge(from_gene, to_gene, weight=weight)

# Define layout (you can try spring_layout or circular_layout)
pos = nx.spring_layout(G, seed=42)

# Get edge weights for thickness
weights = [G[u][v]['weight'] * 5 for u, v in G.edges()]  # scale up for visibility

plt.figure(figsize=(8, 6))
nx.draw(
    G,
    pos,
    with_labels=True,
    node_size=2000,
    node_color="lightgreen",
    edge_color="gray",
    width=weights,
    font_size=12,
    arrows=True,
)
edge_labels = {(u, v): f"{G[u][v]['weight']:.2f}" for u, v in G.edges()}
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="red")

plt.title("Gene Interaction Network (Edge Thickness = Transition Probability)")
plt.show()