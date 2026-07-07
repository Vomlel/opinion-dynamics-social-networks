# Modelování šíření názorů v sociálních sítích



## Abstrakt - CZ

Tato bakalářská práce se zabývá modelováním šíření názorů v sociálních sítích pomocí agentových modelů názorové dynamiky. Sociální síť je reprezentována jako neorientovaný graf, jehož vrcholy představují aktéry a hrany možné interakce mezi nimi. Práce se zaměřuje na vliv parametrů modelů a síťové struktury na vznik konsenzu, polarizace, fragmentace a názorových shluků.

V práci jsou implementovány a analyzovány Deffuantův model a model Hegselmanna–Krauseho. Simulace jsou prováděny na třech typech sítí: náhodné síti Erdős–Rényiho typu, bezškálové síti Barabási–Albertova typu a reálné síti politicky orientovaných facebookových stránek. Výsledky jsou vyhodnocovány pomocí počtu názorových shluků, průměrné párové vzdálenosti názorů a polarizačního indexu Duclos–Esteban–Ray.

Součástí práce je také návrh rozšířeného modelu HRDCR, který doplňuje základní Deffuantův model o strukturální podobnost uzlů, komunitní blízkost a uzlovou setrvačnost. Výsledky ukazují, že mez důvěry významně ovlivňuje přechod mezi fragmentací a konsenzem a že síťová struktura může měnit výsledné rozložení názorů. Rozšířený model HRDCR v experimentech vedl k výraznějšímu lokálnímu shlukování názorů a v doplňkovém porovnání se nejvíce přiblížil empirickému rozložení politických postojů členů 119. Senátu Spojených států.

Klíčová slova: názorová dynamika, sociální sítě, agentové modelování, omezená důvěra, Deffuantův model, Hegselmann–Krause model, polarizace, HRDCR



## Abstract - EN

This bachelor’s thesis focuses on modeling the spread of opinions in social networks using agent-based opinion dynamics models. a social network is represented as an undirected graph whose vertices correspond to actors and whose edges represent possible interactions between them. The thesis examines how model parameters and network structure influence the emergence of consensus, polarization, fragmentation, and opinion clusters.

The thesis implements and analyzes the Deffuant model and the Hegselmann–Krause model. Simulations are performed on three types of networks: an Erdős–Rényi random network, a Barabási–Albert scale-free network, and a real network of politically oriented Facebook pages. The results are evaluated using the number of opinion clusters, the average pairwise opinion distance, and the Duclos–Esteban–Ray polarization index.

The thesis also proposes an extended HRDCR model, which adds structural node similarity, community proximity, and node inertia to the basic Deffuant model. The results show that the confidence threshold significantly affects the transition between fragmentation and consensus and that network structure can influence the final opinion distribution. In the experiments, the extended HRDCR model produced stronger local opinion clustering and, in a supplementary comparison, came closest to the empirical distribution of political attitudes among members of the 119th United States Senate.

Keywords: opinion dynamics, social networks, agent-based modeling, bounded confidence, Deffuant model, Hegselmann–Krause model, polarization, HRDCR


## Reprodukce experimentů

Projekt je implementován v jazyce Julia a byl ověřen s Julií 1.11. Závislosti jsou uvedeny v `code/Project.toml`. Po naklonování repozitáře je lze nainstalovat příkazem:

```bash
julia --project=code -e 'using Pkg; Pkg.instantiate()'
```

Hlavní experiment se z kořene projektu spouští takto:

```bash
julia --project=code code/main.jl -c config.toml
```

Konfigurační cesta je vyhodnocena relativně ke složce `code`, takže uvedený příkaz používá [code/config.toml](code/config.toml). Výchozí konfigurace je produkční a může běžet dlouho. Před demonstračním spuštěním je vhodné snížit `runs`, počet iterací a počet testovaných parametrů nebo použít kopii konfigurace.

Automatické testy lze spustit samostatně:

```bash
julia --project=code code/tests/runtests.jl
```

Analýza empirických hodnot `nominate_dim1` používá data v `data/konkresy/HSall_members.csv`:

```bash
julia --project=code code/analysis/nominate_dim1_analysis.jl
```

## Struktura projektu

```text
code/
  analysis/    samostatné analýzy kongresových dat
  core/        konfigurace a společné pomocné funkce
  metrics/     statické a dynamické charakteristiky
  models/      modely názorové dynamiky
  networks/    generování a načítání sítí
  opinions/    generování počátečních názorů
  reporting/   agregace, grafy a výstupy pro práci
  tests/       automatické testy
  main.jl      hlavní experimentální orchestrátor
data/          vstupní data
out/           generované výsledky (ignorované Gitem)
```

Složka `old_code/` obsahuje historické prototypy a není součástí aktuální implementace. Aktuální kód je výhradně ve složce `code/`.

## Výstupy a reprodukovatelnost

Každý experiment zapisuje výsledky do adresáře určeného položkou `experiment.out_dir` v konfiguraci. Důležité části výstupu jsou:

- `aggregated/static_summary.csv` – agregované vlastnosti sítí,
- `aggregated/dynamic_summary.csv` – agregovaný vývoj dynamických metrik,
- `plots/` – standardní diagnostické grafy,
- `thesis_figures/` – grafy připravené pro sazbu práce,
- `config_used.toml` – kopie konfigurace konkrétního běhu.

Náhodnost je řízena hodnotou `seed_base`; jednotlivé kombinace modelu, sítě a opakování dostávají deterministicky odvozené seedy. Výstupní složka `out/` není verzována, protože produkční běhy mohou obsahovat velké množství dat.

## Použité modely a data

Aktuální experimentální část obsahuje Deffuantův model, model Hegselmann–Krause a navržené strukturální rozšíření HRDCR. Experimenty pracují se sítěmi Erdős–Rényi, Barabási–Albert a s reálnou sítí facebookových stránek politiků. Soubor `HSall_members.csv` slouží k doplňkovému porovnání s empirickou první dimenzí DW-NOMINATE.

Text práce a pracovní podklady jsou ve složce `doc/`. Hlavní zdrojové soubory dat jsou ve složce `data/`.
