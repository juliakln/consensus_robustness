Expected times - use PRISM, in options: change "Engine" to 'sparse', increase value of "Termination max. iterations"; then verify property

Reaching time: 
- reward: !(((x+Zx>=m)&(x>=y+d)) | ((y+Zy>=m)&(y>=x+d))) : 1; --> for not being in a majority state
- property: R=? [ F ("x_wins"|"y_wins") ]; --> for reaching a majority state

Holding time:
- reward: (((x+Zx>=m)&(x>=y+d)) | ((y+Zy>=m)&(y>=x+d))) : 1; --> for being in a majority state
- property:  R=? [ F (("x_wins"|"y_wins")&((!("y_wins")|(F (!("y_wins"))))&(!("x_wins")|(F (!("x_wins")))))) ]; --> for leaving a majority state after reaching it