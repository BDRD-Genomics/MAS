$('button.button').click(function() {
    var load = document.getElementById("load");
    var name = document.forms["Form"]["name"].value;
    var upload = document.forms["Form"]["upload"].value;
    var repeat = document.forms["Form"]["upload"].value;

    if(load.style.display === "none" && name.length !== 0 && repeat.length !== 0 && !upload.isNull) {
        load.style.display = "block";
    }else{
        load.style.display = "none";
    }

});
