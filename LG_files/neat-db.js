var NeatDataBase = function() {
  var dbName = "neat";
  var dbVersion = 1.0;
  var neatDB = {};

  var indexedDB = window.indexedDB || 
                  window.webkitIndexedDB || 
                  window.mozIndexedDB;
  neatDB.indexedDB = {};
  neatDB.indexedDB.db = {};

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
});
