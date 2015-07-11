var AI = function() {
  // nothing
};

AI.prototype.decide = function(game, map, brick, bar) {
  var getColor = function(value) {
    switch (value) {
      case 1: case 3:
        return 1;
      case 2: case 4:
        return 2;
    }
  };
  // we'd like to compute current "goodness"
  // which is determined by each 2 by 2 squares on the map,
  //  4  2  1 -1
  // xx xo xo xo
  // xx xx xo ox
  // @col and @row is the position of bottom left corner of 2 by 2
  // we assume that @col and @row are valid values.
  var computeGoodness = function(col, row) {
    var grid = map.grid;
    var nSame = 0;
    for (var dc = 0; dc <= 1; ++ dc) {
      for (var dr = 0; dr <= 1; ++ dr) {
        if (grid[col + dc][row + dr] == 0) return 0; // has empty tile
        if (getColor(grid[col][row]) == getColor(grid[col + dc][row + dr])) {
          nSame ++;
        }
      }
    }
    switch (nSame) {
      case 4:
        return 4;
      case 1: case 3:
        return 2;
      case 2:
        if (getColor(grid[col][row]) == getColor(grid[col + 1][row + 1])) {
          return -1;
        } else {
          return 1;
        }
    }
    return 0;
  };

  var COL_NUM = Board.COL_NUM;
  var ROW_NUM = Board.ROW_NUM;

  var bestValue = -1e8, maxCol = -1, maxTurn = -1;
  for (var col = 0; col < COL_NUM - 1; ++ col) {
    // if we drop the brick at column "col", then the left half of brick would
    // be placed at
    //    map.colHeight[col] + 2 and map.colHeight[col] + 1
    if (map.colHeight[col] >= ROW_NUM - 1 ||
        map.colHeight[col + 1] >= ROW_NUM - 1) {
      // not a valid option
      continue;
    }
    // make a copy, prevent editing, order: 0 3
    //                                      1 2
    var color = brick.color.slice(0, 4); 
    var hL = map.colHeight[col], hR = map.colHeight[col + 1];
    for (var turn = 0; turn < 4; ++ turn) {
      map.grid[col][hL] = color[1] + 1;
      map.grid[col][hL + 1] = color[0] + 1;
      map.grid[col + 1][hR] = color[2] + 1;
      map.grid[col + 1][hR + 1] = color[3] + 1;

      // +--+
      // |**| function "computeGoodness" would compute the value of 
      // |x*| this 2 by 2, where parameter col, row is position of bottom
      // +--+ left cornor.
      //      now, if we add a brick at position c, r, then, the diffrerence of
      //  ox  total goodness would be goodness of
      //  oo  (c, r), (c - 1, r), (c - 1, r - 1), (c, r - 1)

      var tiles = [[col, hL], [col, hL + 1], [col + 1, hR], [col + 1, hR + 1]];
      var computationRequired = [];
      for (var dc = -1; dc <= 0; ++ dc) {
        for (var dr = -1; dr <= 0; ++ dr) {
          for (var i in tiles) {
            var tile = tiles[i];
            computationRequired[JSON.stringify([tile[0] + dc, tile[1] + dr])] = true;
          }
        }
      }

      var dGoodness = 0;
      for (var i in computationRequired) {
        var tile = JSON.parse(i);
        var c = tile[0];
        var r = tile[1];
        if (c < 0 || c + 1 >= COL_NUM || r < 0 || r + 1 >= ROW_NUM) continue;
        dGoodness += computeGoodness(c, r);
      }

      if (dGoodness > bestValue) {
        bestValue = dGoodness;
        maxCol = col;
        maxTurn = turn;
      }
      
      // recover
      map.grid[col][hL] = 0;
      map.grid[col][hL + 1] = 0;
      map.grid[col + 1][hR] = 0;
      map.grid[col + 1][hR + 1] = 0;

      // simulate pressing 'up' key
      var tmp = color.slice(1, 4);
      tmp.push(color[0]);
      color = tmp;
    }
  }
  if (maxCol == -1) {
    return -1;
  }
  if (maxCol < brick.col) {
    return 'left';
  }
  if (maxCol > brick.col) {
    return 'right';
  }
  if (maxTurn != 0) {
    return 'up';
  }
  return 'down';
};

var computerAI = new AI();

$(document).ready(function() {
  $('#AI_button').click(function() {
    if (computerAI.enabled) {
      computerAI.enabled = false;
      $('#AI_button').text('Start Demo');
    } else {
      computerAI.enabled = true;
      $('#AI_button').text('Stop Demo');
    }
  });
});
