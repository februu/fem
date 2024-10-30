```cpp
dN_dKsi[0][0] = -0.25 * (1 - gaussNodes[1][0]);
dN_dKsi[0][1] = 0.25 * (1 - gaussNodes[1][0]);
dN_dKsi[0][2] = 0.25 * (1 + gaussNodes[1][0]);
dN_dKsi[0][3] = -0.25 * (1 + gaussNodes[1][0]);

dN_dKsi[1][0] = -0.25 * (1 - gaussNodes[1][0]);
dN_dKsi[1][1] = 0.25 * (1 - gaussNodes[1][0]);
dN_dKsi[1][2] = 0.25 * (1 + gaussNodes[1][0]);
dN_dKsi[1][3] = -0.25 * (1 + gaussNodes[1][0]);

dN_dKsi[2][0] = -0.25 * (1 - gaussNodes[1][1]);
dN_dKsi[2][1] = 0.25 * (1 - gaussNodes[1][1]);
dN_dKsi[2][2] = 0.25 * (1 + gaussNodes[1][1]);
dN_dKsi[2][3] = -0.25 * (1 + gaussNodes[1][1]);

dN_dKsi[3][0] = -0.25 * (1 - gaussNodes[1][1]);
dN_dKsi[3][1] = 0.25 * (1 - gaussNodes[1][1]);
dN_dKsi[3][2] = 0.25 * (1 + gaussNodes[1][1]);
dN_dKsi[3][3] = -0.25 * (1 + gaussNodes[1][1]);

// Eta
dN_dEta[0][0] = -0.25 * (1 - gaussNodes[1][0]);
dN_dEta[0][1] = -0.25 * (1 + gaussNodes[1][0]);
dN_dEta[0][2] = 0.25 * (1 + gaussNodes[1][0]);
dN_dEta[0][3] = 0.25 * (1 - gaussNodes[1][0]);

dN_dEta[1][0] = -0.25 * (1 - gaussNodes[1][1]);
dN_dEta[1][1] = -0.25 * (1 + gaussNodes[1][1]);
dN_dEta[1][2] = 0.25 * (1 + gaussNodes[1][1]);
dN_dEta[1][3] = 0.25 * (1 - gaussNodes[1][1]);

dN_dEta[2][0] = -0.25 * (1 - gaussNodes[1][1]);
dN_dEta[2][1] = -0.25 * (1 + gaussNodes[1][1]);
dN_dEta[2][2] = 0.25 * (1 + gaussNodes[1][1]);
dN_dEta[2][3] = 0.25 * (1 - gaussNodes[1][1]);

dN_dEta[3][0] = -0.25 * (1 - gaussNodes[1][0]);
dN_dEta[3][1] = -0.25 * (1 + gaussNodes[1][0]);
dN_dEta[3][2] = 0.25 * (1 + gaussNodes[1][0]);
dN_dEta[3][3] = 0.25 * (1 - gaussNodes[1][0]);
```
