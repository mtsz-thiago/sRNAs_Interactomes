import neo4j
import networkx as nx
import glob
from tqdm import tqdm


uri = "bolt://neo4j:7687"
user = "neo4j"
password = "Password"

# Get a list of all gml files in the directory
gml_files = glob.glob("../output/interactome_graphs/*.gml")

# Function to load a gml file into Neo4j
def load_gml_to_neo4j(file, driver):
    # Load the gml file into a NetworkX graph
    graph = nx.read_gml(file)

    # Convert the graph to a Neo4j graph
    with driver.session() as session:
        for node in graph.nodes(data=True):
            session.run("CREATE (n:Node {id: $id, attr: $attr})", id=node[0], attr=node[1])

        for edge in graph.edges(data=True):
            session.run("""
                MATCH (a:Node {id: $source})
                MATCH (b:Node {id: $target})
                CREATE (a)-[:CONNECTED_TO {attr: $attr}]->(b)
            """, source=edge[0], target=edge[1], attr=edge[2])
            
if __name__ == "__main__":
  
    with neo4j.GraphDatabase.driver(uri=uri, auth=(user, password)) as driver:
        
        # Verify that we can connect to the database
        driver.verify_connectivity()
        print('driver ok')
        
        for file in tqdm(gml_files):    
            load_gml_to_neo4j(file, driver)
            print(f"Loaded {file} into Neo4j")