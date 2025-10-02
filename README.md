# BioCompare

> Quickly explore and compare proteins from multiple trusted sources at one place!

This is a full-stack bioinformatics web application designed for the rapid comparison of two proteins using their Ensembl IDs. It aggregates data from major biological databases and presents the comparison through interactive visualizations and tables.

## 🚀 Live Demo

(https://biocompare-by-rudra.onrender.com/)

## ✨ Key Features

* **Multi-Source Data Aggregation**: The application fetches and integrates data from three primary sources: Ensembl, UniProt, and NCBI Gene.
* **Interactive Venn Diagrams**: Compares Gene Ontology (GO) IDs and GO Terms between two proteins and visualizes the intersection and differences.
* **Detailed Gene Identity View**: Displays a summary card for each protein, showing its gene symbol, full name, and organism.
* **Tabular Data View**: Presents the raw, aggregated GO data in a clear table format.
* **CSV Export**: Allows users to download the comparison data as a CSV file.
* **Robust Backend**: The Flask backend features a resilient API client with retry and backoff logic to handle network errors or API rate limits.
* **Responsive Design**: The user interface is designed to be functional across different screen sizes.

## 🛠️ Technology Stack

* **Backend**: Python, Flask, Flask-CORS
* **Frontend**: HTML5, CSS3, vanilla JavaScript
* **Data Visualization**: d3.js, venn.js
* **API Interaction**: `requests` library in Python
* **Deployment**: Monolithic Flask app deployed on Render

## ⚙️ How It Works: The Data Pipeline

1.  **Input**: The user provides two Ensembl IDs via the web interface.
2.  **Primary Fetch**: The Flask backend queries the Ensembl REST API for primary gene info and GO terms.
3.  **Data Enrichment**: The backend then uses fallback mechanisms to query UniProt and NCBI to find corresponding entries and gather additional GO terms.
4.  **Data Merging**: GO annotations from all sources are merged into a single, comprehensive list for each protein.
5.  **Response**: The final, aggregated JSON is sent to the frontend for visualization.

## Usage

1.  Navigate to the live demo URL.
2.  Enter two Ensembl IDs (or use the "Copy example" button to see a proper example).
3.  Click "Submit".
4.  Use the controls on the left to switch between the **GO IDs**, **GO Terms**, **Gene Info**, and **Table** comparison views.
