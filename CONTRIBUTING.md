# Contributing Guide - UVR Hybrid Modelling v2 - Skvader-1 Flight Hybrid

This repository is a shared development fork of the **King_Hybrid** simulation framework.  
To avoid merge conflicts and broken commits, please follow these rules carefully.

---

## 🔒 Branch Model

| Branch | Purpose | Who Can Push |
|---------|----------|--------------|
| `main` | Stable, tested release code | **Protected** — only maintainers merge |
| `dev` | Integration / staging branch | All contributors (via Pull Requests) |
| `feature/<name>` | Individual development branches | Each contributor |

---

## 🧩 Workflow Overview

1. **Fork this repo** into your own GitHub account.  
   - Example: `github.com/teammate/uvr_hybrid_modelling_v2`

2. **Create a new branch** for your work:
   ```bash
   git checkout -b feature/<your-name>

3.**Make a new Skvader-1 input file**
    Click --> Skvader_1_Flight_Hybrid
    ctrl+c, ctrl+v
    Right click the copy and rename:
    Skvader_1_Flight_Hybrid_[your name]

    copy/paste from other files if you want to keep your file in sync with others
    this is the easiest way i can think of

### Happy Modelling and Engine Design!