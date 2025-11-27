import csv
import sys
import random

# Importa o algoritmo do outro arquivo
try:
    from algoritmo_dijkstra import dijkstra_sinalizacao
except ImportError:
    print("ERRO: Não foi possível importar 'algoritmo_dijkstra.py'.")
    print("Verifique se os dois arquivos estão na mesma pasta.")
    sys.exit()

def carregar_dados_string(caminho_arquivo, limite_confianca):
    """
    Lê o arquivo do STRING DB e monta o grafo ponderado.
    """
    grafo = {}
    print(f"Lendo base de dados: {caminho_arquivo} (com limite de confiança {limite_confianca})")
    
    try:
        with open(caminho_arquivo, 'r') as f:
            leitor = csv.reader(f, delimiter=' ')
            next(leitor) # Pula cabeçalho
            
            contagem = 0
            for linha in leitor:
                # Tratamento para linhas incompletas
                if len(linha) < 3: continue
                
                p1, p2 = linha[0], linha[1]
                try:
                    score = float(linha[2])
                except ValueError:
                    continue # Pula se o score não for número
                
                # Filtra interações fracas
                if score < limite_confianca: continue
                
                # Conversão Score -> Peso (Custo)
                # Quanto maior o score, menor o custo (caminho mais curto/fácil)
                peso = 1.0 - (score / 1000.0)
                
                if p1 not in grafo: grafo[p1] = {}
                if p2 not in grafo: grafo[p2] = {}
                
                # Armazena a aresta (grafo não direcionado)
                if p2 not in grafo[p1] or peso < grafo[p1][p2]:
                    grafo[p1][p2] = peso
                    grafo[p2][p1] = peso
                    contagem += 1
                    
        print(f"\nSucesso! Grafo montado com {len(grafo)} proteínas e {contagem} conexões.")
        return grafo

    except FileNotFoundError:
        print(f"ERRO: O arquivo '{caminho_arquivo}' não foi encontrado na pasta.")
        return {}
    
def encontrar_par(grafo):
    """
    Escolhe uma origem e um destino aleatórios que tenham um caminho entre si.
    1. Pega uma origem aleatória.
    2. Vê quem ela alcança (BFS rápido).
    3. Escolhe um destino aleatório dentre os alcançáveis.
    """
    print("\nSorteando um par com conexão válida...")
    
    lista_proteinas = list(grafo.keys())
    random.shuffle(lista_proteinas) # Embaralha para ser sempre diferente
    
    for origem in lista_proteinas:
        # Faz uma busca rápida (BFS) para ver quem essa origem alcança
        fila = [origem]
        visitados = {origem}
        alcançaveis_distantes = [] # Lista de candidatos a destino
        
        # Mapa de distâncias
        distancia = {origem: 0}
        
        idx = 0
        while idx < len(fila):
            atual = fila[idx]
            idx += 1
            
            dist_atual = distancia[atual]
            
            # Se a distância for razoável (ex: >= 2 passos), é um bom candidato a destino
            if dist_atual >= 2:
                alcançaveis_distantes.append(atual)
            
            # Limita a busca para não travar em grafos gigantes
            if len(alcançaveis_distantes) > 20: 
                break 
            
            for vizinho in grafo[atual]:
                if vizinho not in visitados:
                    visitados.add(vizinho)
                    distancia[vizinho] = dist_atual + 1
                    fila.append(vizinho)
        
        # Se essa origem consegue chegar em lugares distantes, escolhe um destino aleatório
        if alcançaveis_distantes:
            destino = random.choice(alcançaveis_distantes)
            print(f"\nPar encontrado após busca: {origem} -> {destino}")
            return origem, destino

    # Se não achar nenhum par, escolhe aleatoriamente
    p1 = random.choice(lista_proteinas)
    p2 = random.choice(lista_proteinas)
    print(f"\nNenhum par com conexão encontrado. Usando par aleatório: {p1} -> {p2}")
    return p1, p2

# Lógica principal
if __name__ == "__main__":

    # Nome do arquivo baixado
    arquivo_dados = "511145.protein.links.v12.0.txt" 
    
    # Limite 700 considera apenas "Alta Confiança". Se não achar caminhos, baixar pra 400.
    limite = 700 
    
    # Carrega os dados
    grafo_real = carregar_dados_string(arquivo_dados, limite_confianca=limite)
    
    if grafo_real:
        # Escolha aleatória de origem e destino 
        origem, destino = encontrar_par(grafo_real)
        
        print(f"\n--- Calculando Via de Sinalização ---")
        print(f"Origem (Receptor): {origem}")
        print(f"Destino (Alvo):    {destino}")
        
        # Executa o algoritmo de Dijkstra
        caminho, custo = dijkstra_sinalizacao(grafo_real, origem, destino)
        
        if caminho:
            print("\n" + "="*60)
            print(f"VIA DE SINALIZAÇÃO ENCONTRADA")
            print("="*60)
            print(f"Passos (Nós): {len(caminho)}")
            print(f"Custo total (Resistência): {custo:.4f}")
            print("-" * 60)
            print(" -> ".join(caminho))
            print("-" * 60)
            print("\n>>> PARA O CYTOSCAPE:")
            print(" OR ".join(caminho))
            print("-" * 60)
        else:
            print("\nNão foi possível encontrar um caminho (tente rodar novamente, ou diminuir o limite de confiança).")