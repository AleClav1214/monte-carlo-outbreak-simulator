"""Transmission network visualization."""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx

from outbreak_simulator.transmission_inference.epi_reconstruction import ReconstructedTransmissionNetwork


def plot_transmission_network(
    network: ReconstructedTransmissionNetwork,
    title: str = "Reconstructed transmission network",
    highlight_top_n_spreaders: int = 3,
) -> plt.Figure:
    """Draw the reconstructed who-infected-whom network, with edge width
    proportional to reconstruction confidence (P(i infected j)) and node
    size/color proportional to each case's reconstructed onward-transmission
    (superspreader) score."""
    graph = network.graph
    fig, ax = plt.subplots(figsize=(9, 7))

    if graph.number_of_nodes() == 0:
        ax.text(0.5, 0.5, "No cases to display", ha="center", va="center")
        return fig

    pos = nx.spring_layout(graph, seed=42, k=1.5 / max(graph.number_of_nodes() ** 0.5, 1))
    scores = network.superspreader_scores
    top_spreaders = set(n for n, _ in sorted(scores.items(), key=lambda x: -x[1])[:highlight_top_n_spreaders])

    node_sizes = [200 + 400 * scores.get(n, 0) for n in graph.nodes]
    node_colors = ["#c0392b" if n in top_spreaders else "#2980b9" for n in graph.nodes]
    edge_widths = [1 + 3 * graph[u][v]["weight"] for u, v in graph.edges]
    edge_alphas = [max(0.15, min(graph[u][v]["weight"], 1.0)) for u, v in graph.edges]

    nx.draw_networkx_nodes(graph, pos, node_size=node_sizes, node_color=node_colors, alpha=0.85, ax=ax)
    for (u, v), width, alpha in zip(graph.edges, edge_widths, edge_alphas):
        nx.draw_networkx_edges(graph, pos, edgelist=[(u, v)], width=width, alpha=alpha, edge_color="#7f8c8d",
                                arrows=True, arrowsize=10, connectionstyle="arc3,rad=0.05", ax=ax)
    if graph.number_of_nodes() <= 60:
        nx.draw_networkx_labels(graph, pos, font_size=7, ax=ax)

    ax.set_title(f"{title}\n(red = top {highlight_top_n_spreaders} reconstructed superspreaders, edge width = confidence)")
    ax.axis("off")
    fig.tight_layout()
    return fig


def plot_generation_size_bar(generation_sizes, title: str = "Cases per generation") -> plt.Figure:
    """Simple bar chart of the branching-process model's generation_sizes output."""
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(range(len(generation_sizes)), generation_sizes, color="#16a085")
    ax.set_xlabel("Generation")
    ax.set_ylabel("New cases")
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig
