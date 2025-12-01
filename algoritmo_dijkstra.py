import heapq

def dijkstra_sinalizacao(grafo, proteina_origem, proteina_destino):
    """
    Implementa o algoritmo de Dijkstra para encontrar a via de sinalização
    (caminho de menor custo) entre duas proteínas.
    """
    
    # Inicialização
    # Distâncias começam como infinito, exceto a origem que é 0
    distancias = {nodo: float('infinity') for nodo in grafo}
    distancias[proteina_origem] = 0
    
    # Dicionário para reconstruir o caminho (quem veio antes de quem)
    predecessores = {nodo: None for nodo in grafo}
    
    # Fila de prioridade: armazena tuplas (custo_acumulado, proteina_atual)
    # O heapq sempre coloca o menor custo no topo
    fila_prioridade = [(0, proteina_origem)]
    
    while fila_prioridade:
        # Remove a proteína com o menor custo atual da fila
        custo_atual, nodo_atual = heapq.heappop(fila_prioridade)
        
        # Se chegou ao destino, para o loop
        if nodo_atual == proteina_destino:
            break
        
        # Se o custo retirado da fila for maior que o já conhecido, ignora
        if custo_atual > distancias[nodo_atual]:
            continue
        
        # Verifica todos os vizinhos (proteínas que interagem)
        for vizinho, peso_aresta in grafo[nodo_atual].items():
            distancia_nova = custo_atual + peso_aresta
            
            # Se encontrou um caminho mais curto para o vizinho, atualiza
            if distancia_nova < distancias[vizinho]:
                distancias[vizinho] = distancia_nova
                predecessores[vizinho] = nodo_atual
                heapq.heappush(fila_prioridade, (distancia_nova, vizinho))
    
    # Reconstrução do caminho (Backtracking)
    caminho = []
    nodo_atual = proteina_destino
    
    # Se o destino continua com distância infinita, não há caminho
    if distancias[proteina_destino] == float('infinity'):
        return None, float('infinity')
        
    while nodo_atual is not None:
        caminho.insert(0, nodo_atual)
        nodo_atual = predecessores[nodo_atual]
        
    return caminho, distancias[proteina_destino]

# --- EXEMPLO DE EXECUÇÃO ---

if __name__ == "__main__":
    # Simulação de uma Rede de Interação Proteína-Proteína (PPI).
    # IMPORTANTE: No Dijkstra, queremos o MENOR custo. 
    # Em biologia, interações têm "scores" de confiança (ex: 0.9 é forte).
    # Para usar Dijkstra, transformamos o score em custo/peso: 
    # Peso = 1.0 - Score (ou 1/Score). Assim, maior confiança = menor custo.
    
    grafo_interacoes = {
        'Proteina_A': {'Proteina_B': 0.1, 'Proteina_C': 0.5}, # 0.1 = alta afinidade (score 0.9)
        'Proteina_B': {'Proteina_A': 0.1, 'Proteina_D': 0.2, 'Proteina_E': 0.8},
        'Proteina_C': {'Proteina_A': 0.5, 'Proteina_F': 0.4},
        'Proteina_D': {'Proteina_B': 0.2, 'Proteina_G': 0.1}, # Caminho forte
        'Proteina_E': {'Proteina_B': 0.8, 'Proteina_G': 0.9}, # Caminho fraco
        'Proteina_F': {'Proteina_C': 0.4, 'Proteina_G': 0.6},
        'Proteina_G': {'Proteina_D': 0.1, 'Proteina_E': 0.9, 'Proteina_F': 0.6}
    }

    inicio = 'Proteina_A'  # Ex: Receptor de membrana
    fim = 'Proteina_G'     # Ex: Fator de transcrição no núcleo

    print(f"--- Buscando Via de Sinalização: {inicio} -> {fim} ---")
    
    caminho, custo = dijkstra_sinalizacao(grafo_interacoes, inicio, fim)

    if caminho:
        print(f"Via encontrada: {' -> '.join(caminho)}")
        print(f"Custo total do caminho (resistência): {custo}")
        print("Interpretação: Esta é a cascata de sinalização mais provável.")
    else:
        print("Não foi encontrada nenhuma via de interação entre as proteínas.")