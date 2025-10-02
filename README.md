BioCompare
BioCompare is a full-stack bioinformatics web application designed for the rapid comparison of two proteins using their Ensembl IDs. It aggregates data from major biological databases and presents the comparison through interactive visualizations and tables.

🚀 Live Demo
(https://biocompare-by-rudra.onrender.com/)

✨ Key Features
Multi-Source Data Aggregation: The application fetches and integrates data from three primary sources: Ensembl, UniProt, and NCBI Gene.

Interactive Venn Diagrams: Compares Gene Ontology (GO) IDs and GO Terms between two proteins and visualizes the intersection and differences with Venn diagrams powered by d3.js and venn.js.

Detailed Gene Identity View: Displays a summary card for each protein, showing its gene symbol, full name, organism, and Ensembl ID.

Tabular Data View: Presents the raw, aggregated GO data in a clear table format.

CSV Export: Allows users to download the comparison data from the table view as a CSV file.

Robust Backend: The Flask backend features a resilient API client with retry and backoff logic to handle network errors or API rate limits gracefully.

Responsive Design: The user interface is designed to be accessible and functional across different screen sizes, including mobile devices.

🛠️ Technology Stack
Backend: Python, Flask, Flask-CORS.

Frontend: HTML5, CSS3, vanilla JavaScript.

Data Visualization: d3.js, venn.js.

API Interaction: The requests library in Python is used for server-side API calls.

Deployment: The application is a monolithic Flask app that serves both the API and the static frontend.

⚙️ How It Works: The Data Pipeline
The application's backend is the core of its functionality, performing a multi-stage data aggregation process for each user request:

Input: The user provides two Ensembl IDs via the web interface.

Primary Fetch: The Flask backend receives the IDs and first queries the Ensembl REST API to get primary gene information (symbol, description, organism) and a list of GO term cross-references.

UniProt Fallback: To enrich the data, the backend uses the Ensembl ID to find a corresponding UniProt entry. If found, it fetches additional GO terms from the UniProt REST API.

NCBI Fallback: The backend also queries NCBI E-utils to find a corresponding NCBI Gene ID, further confirming the gene's identity.

Data Merging: GO annotations from Ensembl and UniProt are merged into a single, comprehensive list for each protein, removing duplicates.

Response: The final, aggregated JSON object is sent to the frontend, which then renders the various comparison views (Venn diagrams, tables, etc.).

Usage-
Navigate to the live demo URL.

Enter an Ensembl ID into the "Ensembl ID A" field and another into the "Ensembl ID B" field. You can use the provided "Copy example" button to populate the fields with sample data.

Click the "Submit" button.

After the data is fetched, use the controls on the left to switch between the GO IDs, GO Terms, Gene Info, and Table comparison views.







