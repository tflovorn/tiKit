# coding=UTF-8

def __setup():
    symtext = """α &alpha;   \\alpha     {/Symbol a}    224
    β &beta;    \\beta      {/Symbol b}    225
    &chi &chi;     \\chi       {/Symbol c}
    &Chi &Chi;     \\Chi       {/Symbol C}
    δ &delta;   \\delta     {/Symbol d}    235
    Δ &Delta;   \\Delta     {/Symbol D} 
    ε &epsilon; \\epsilon   {/Symbol e}    238
    Ε &Epsilon; \\Epsilon   {/Symbol E}
    φ &phi;     \\phi       {/Symbol f}    237
    Φ &Phi;     \\Phi       {/Symbol F}    232
    γ &gamma;   \\gamma     {/Symbol g}  
    Γ &Gamma;   \\Gamma     {/Symbol G}    226
    η &eta;     \\eta       {/Symbol h}
    Η &Eta;     \\Eta       {/Symbol H}
    ι &iota;    \\iota      {/Symbol i}
    Ι &Iota;    \\Iota      {/Symbol I}
    κ &kappa;   \\kappa     {/Symbol k}
    Κ &Kappa;   \\Kappa     {/Symbol K}
    &lambda &lambda;  \\lambda    {/Symbol l}
    &Lambda &Lambda;  \\Lambda    {/Symbol L}
    μ &mu;      \\mu        {/Symbol m}    230
    Μ &Mu;      \\Mu        {/Symbol M}
    ν &nu;      \\nu        {/Symbol n}
    Ν &Nu;      \\Nu        {/Symbol N}
    &omicron &omicron; \\omicron   {/Symbol o}
    &Omicron &Omicron; \\Omicron   {/Symbol O}
    π &pi;      \\pi        {/Symbol p}    227
    Π &Pi;      \\Pi        {/Symbol P}
    θ &theta;   \\theta     {/Symbol q}  
    Θ &Theta;   \\Theta     {/Symbol Q}    233
    ρ &rho;     \\rho       {/Symbol r}
    Ρ &Rho;     \\Rho       {/Symbol R}
    σ &sigma;   \\sigma     {/Symbol s}    229
    Σ &Sigma;   \\Sigma     {/Symbol S}    228
    τ &tau;     \\tau       {/Symbol t}    231
    Τ &Tau;     \\Tau       {/Symbol T}
    υ &upsilon; \\upsilon   {/Symbol u}
    Υ &Upsilon; \\Upsilon   {/Symbol U}
    ω &omega;   \\omega     {/Symbol w}
    Ω &Omega;   \\Omega     {/Symbol W}
    ξ &xi;      \\xi        {/Symbol x}
    Ξ &Xi;      \\Xi        {/Symbol X}
    ψ &psi;     \\psi       {/Symbol y}
    Ψ &Psi;     \\Psi       {/Symbol Y}
    ζ &zeta;    \\zeta      {/Symbol z}
    Ζ &Zeta;    \\Zeta      {/Symbol Z}""".splitlines()

    for line in symtext:
        ls = line.split()
        name = ls[1][1:-1]
        symb = ls[3] + " " + ls[4]
        symbolsDict[name] = symb

symbolsDict = {}
__setup()

def Convert(name):
    if name in symbolsDict:
        return symbolsDict[name]
    else:
        return None
