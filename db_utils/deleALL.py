from neo4j import GraphDatabase

class Neo4jService(object):
    def __init__(self, uri, user, password):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self._driver.close()

    def delete_all(self):
        with self._driver.session() as session:
            return session.write_transaction(self._delete_all)

    @staticmethod
    def _delete_all(tx):
        query = """
        MATCH (n)
        DETACH DELETE n
        """
        return tx.run(query)

# Replace 'bolt://localhost:7687', 'neo4j', 'password' with your actual Neo4j credentials
service = Neo4jService('bolt://localhost:7687', 'neo4j', 'Password')
service.delete_all()
service.close()