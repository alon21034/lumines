var NeatDataBase = function() {
  var dbName = "neat";
  var dbVersion = 1.0;
  var neatDB = {};

  var indexedDB = window.indexedDB || 
                  window.webkitIndexedDB || 
                  window.mozIndexedDB;
  neatDB.indexedDB = {};
  neatDB.indexedDB.db = {};

  this.indexedDB = neatDB.indexedDB;
  if ('webkitIndexedDB' in window) {
    window.IDBTransaction = window.webkitIDBTransaction;
    window.IDBKeyRange = window.webkitIDBKeyRange;
  }

  this.open = function() {
    var request = indexedDB.open(dbName, dbVersion);

    request.onsuccess = function(e) {
      console.log ("success to open DB: " + dbName);
      neatDB.indexedDB.db = e.target.result;
      var db = neatDB.indexedDB.db;
    };

    request.onupgradeneeded = function(e) {
      console.log ("upgrade required");
      neatDB.indexedDB.db = e.target.result;
      var db = neatDB.indexedDB.db;
      if (db.objectStoreNames.contains("generation")) {
        db.deleteObjectStore("generation");
      }

      var store = db.createObjectStore("generation",
          {keyPath: "id", autoIncrement: true});

      // store.createIndex("parent_id", "parent_id", {unique: false});
      store.createIndex("generation", "generation", {unique: false});
      store.createIndex("maxFitness", "Fitness", {unique: false});
    };
  };

  this.resetDataBase = function() {
    var transaction = neatDB.indexedDB.db.transaction(["generation"], "readwrite");
    var store = transaction.objectStore("generation");
    var request = store.clear();
    request.onsuccess = function(e) {
      console.log("done!");
    }
  }

  this.saveGeneration = function (pool) {
    // the @pool represents a generation
    var transaction = neatDB.indexedDB.db.transaction(["generation"], "readwrite");

    // transaction.oncomplete = function(e) { };
    transaction.onerror = function(e) {
      console.log(e);
    };

    var store = transaction.objectStore("generation");
    var request = store.add(pool);
    // request.onsuccess = function (e) { };
  };

  this.getGeneration = function(id, callback) {
    var transaction = neatDB.indexedDB.db.transaction(["generation"]);
    var store = transaction.objectStore("generation");
    var request = store.get(id);
    request.onsuccess = function(e) {
      callback(request.result);
    };
    request.onerror = function(e) {
      console.log(e);
    };
  };
};

neatDataBase = new NeatDataBase();
$(document).ready(function() {
  neatDataBase.open();
  $('#update-download-list').click(function() {
    var store = neatDataBase.indexedDB.db.transaction(["generation"]).objectStore("generation");
    $('#download-list').html('');
    store.openCursor().onsuccess = function (e) {
      var cursor = e.target.result;
      if (cursor) {
        $('#download-list').append(
          $("<option></option>").attr("value",cursor.key).text(cursor.key)); 
        cursor.continue();
      }
    };
  });
  var downloadURL;
  $('#create-link').click(function(e) {
    var id;
    $('#download-list option:selected').each(function() {
      id = parseInt($(this).text());
    });
    var link = $('#download-link');
    if (link.attr('href') !== undefined) {
      window.URL.revokeObjectURL(link.attr('href'));
    }

    neatDataBase.getGeneration(id, function(result) {
      var blob = new Blob([result.packed], {type: "octet/stream"});
      link.attr('href', window.URL.createObjectURL(blob));
      link.attr('download', id + ".txt");
      link.text(id);
    });
  });
});

