# lab2ldr

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lkhenayfis.github.io/lab2ldr.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lkhenayfis.github.io/lab2ldr.jl/dev/)
[![Build Status](https://github.com/lkhenayfis/lab2ldr.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lkhenayfis/lab2ldr.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/lkhenayfis/lab2ldr.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/lkhenayfis/lab2ldr.jl)

`lab2ldr` e uma ponte entre o [SDDPlab](https://github.com/rjmalves/sddp-lab) e o
[LinearDecisionRules](https://github.com/bfpc/LinearDecisionRules.jl): a partir de um estudo no
formato do SDDPlab, monta um `LDRModel` (regras de decisao lineares) equivalente ao problema
multiestagio. Com isso, o **mesmo estudo** pode ser resolvido por dois motores — SDDP (via
SDDPlab) e LDR (via lab2ldr) — o que torna o pacote conveniente para ser importado em outros
repositorios (por exemplo, o proprio [CasosLDR](https://github.com/lkhenayfis/CasosLDR)) e
comparar as duas abordagens sobre um mesmo caso.

O repositorio traz dois estudos de exemplo em `data/`:

- `data/hydrotoy2` — problema hidrotermico de 2 estagios, 1 usina hidreletrica;
- `data/hydrotoy3` — o mesmo problema com 3 estagios.

> **Nota:** atualmente apenas o caminho com uma unica usina hidreletrica esta funcional.

## Instalacao

`lab2ldr` nao esta registrado no General e depende de duas bibliotecas em versoes de
desenvolvimento (branches git) que tambem nao estao registradas. Por isso, adicione primeiro
essas dependencias nas revisoes corretas e so entao o `lab2ldr`:

```julia
using Pkg
Pkg.add(url = "https://github.com/rjmalves/sddp-lab.git", rev = "abstract-engine")
Pkg.add(url = "https://github.com/bfpc/LinearDecisionRules.jl.git", rev = "Multistage")
Pkg.add(url = "https://github.com/lkhenayfis/lab2ldr.jl", rev = "v1.0.0")
```

A ordem importa: como essas revisoes nao estao em nenhum registro, precisam ja estar presentes no
ambiente antes da resolucao do `lab2ldr` — do contrario a resolucao tentaria buscar versoes
registradas incompativeis (`SDDPlab` 2.0.0 e a `LinearDecisionRules` do branch `Multistage`). Para
desenvolvimento local, use `Pkg.develop(path = "...")` no lugar de `Pkg.add`.

## Preparando o estudo: chaves `cuts` e `truncation`

O formato de estudo do SDDPlab ainda nao modela nativamente a funcao de custo futuro (condicao de
contorno terminal) nem os limites de truncamento das afluencias. Como o motor LDR precisa dos
dois, o `lab2ldr` os le de duas chaves adicionais no `main.jsonc`, no mesmo nivel de `inputs` e
`engine`:

- **`cuts`** — os cortes lineares que descrevem o custo futuro sobre o volume armazenado no ultimo
  estagio. Sao um `intercept` (vetor, um valor por corte) e `coefficients` (uma linha por corte,
  uma coluna por usina); cada corte `k` impoe `theta >= intercept[k] + coefficients[k, :] .
volume_final`. Sem esse termo terminal o problema fica sem valor da agua e a solucao e
  degenerada (esvaziar o reservatorio nao tem custo).

- **`truncation`** — os limites `[inferior, superior]` de cada afluencia incerta, por usina (`id`)
  e por estagio (`stage`). O conjunto `Uncertainty` do LinearDecisionRules exige suporte finito,
  entao as marginais dos estagios estocasticos (do 2o em diante) sao truncadas nesses limites
  antes de entrarem no modelo. E preciso uma entrada de `bounds` para cada estagio estocastico.

```jsonc
"cuts": {
    "coefficients": [
        [-300.0, -600.0, -1200.0, -2500.0, -7800.0]
    ],
    "intercept": [75000.0, 135000.0, 210000.0, 310000.0, 390000.0]
},
"truncation": [
    {
        "id": 1,
        "bounds": [
            { "stage": 2, "value": [0.0, 120.0] },
            { "stage": 3, "value": [0.0, 120.0] }
        ]
    }
]
```

> **Importante — afluencia do 1o estagio.** O primeiro estagio e deterministico. O `lab2ldr` fixa
> a afluencia do estagio 1 na **media** da marginal daquele estagio (que **nao** e truncada). Para
> garantir que os dois motores resolvam exatamente o mesmo problema, defina a marginal do estagio
> 1 (`season` 1), para toda usina, como uma `Normal` com **media igual ao valor desejado e
> variancia zero** — no arquivo de cenarios de afluencia, `"kind": "Normal", "parameters": [valor,
0]`. Assim o LDR usa a media (o proprio valor) e o SDDP, ao amostrar, so pode obter esse mesmo
> valor: os dois enxergam a mesma condicao inicial deterministica. Uma variancia positiva no
> estagio 1 tornaria o primeiro estagio estocastico apenas para o SDDP, divergindo do LDR.

## Motor LDR (lab2ldr)

`build_model` le o estudo e devolve um `LDRModel` pronto para otimizar. `pkgdir` e usado apenas
para localizar os estudos de exemplo distribuidos com o pacote — em um caso proprio, basta passar
o caminho do diretorio do estudo.

```julia
using lab2ldr
using LinearDecisionRules
using HiGHS

for estudo in ["hydrotoy2", "hydrotoy3"]
    path = pkgdir(lab2ldr, "data", estudo)

    # n_breakpoints controla a aproximacao por partes das regras de decisao das afluencias
    m = build_model(path, HiGHS.Optimizer; n_breakpoints = 3)
    LinearDecisionRules.optimize!(m)

    primal = LinearDecisionRules.objective_value(m)               # limite superior (regra primal)
    dual   = LinearDecisionRules.objective_value(m; dual = true)  # limite inferior (regra dual)
    println("$estudo: primal = $primal | dual = $dual")
end
```

Para ler os coeficientes da regra de decisao de uma variavel, use `get_decision`. Os nomes das
variaveis do sistema sao constantes do SDDPlab (submodulo `Lab`):

```julia
using lab2ldr, LinearDecisionRules, HiGHS
import SDDPlab.Lab: STORED_VOLUME

m = build_model(pkgdir(lab2ldr, "data", "hydrotoy2"), HiGHS.Optimizer)
LinearDecisionRules.optimize!(m)

# termo constante da regra do volume armazenado no 1o estagio
v1 = LinearDecisionRules.get_decision(m, m[STORED_VOLUME][1, 1])
```

## Motor SDDP (SDDPlab), mesmo estudo

O mesmo diretorio de estudo e resolvido pelo motor SDDP do proprio SDDPlab:

```julia
using lab2ldr   # apenas para localizar os estudos de exemplo via pkgdir
using SDDPlab
using HiGHS

for estudo in ["hydrotoy2", "hydrotoy3"]
    path = pkgdir(lab2ldr, "data", estudo)

    study  = SDDPlab.read_study(path)
    model  = SDDPlab.build(study, HiGHS.Optimizer)
    policy = SDDPlab.train(study, model)       # treina a politica SDDP
    sim    = SDDPlab.simulate(study, model)    # simula a politica treinada
end
```
