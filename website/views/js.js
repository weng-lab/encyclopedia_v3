function activeAssembly(){
    var assembly = "hg19";
    if($("#assembly_mm9").is(":checked")){
        assembly = "mm9";
    } else if($("#assembly_mm10").is(":checked")){
        assembly = "mm10";
    }
    return assembly;
}

function updateTissues(){
    $(".selectTissues").hide();

    var assembly = activeAssembly()

    var assays = "H3K27ac";
    if($("#assaysBoth").is(":checked")){
        assays = "Both";
    } else if($("#assaysDNase").is(":checked")){
        assays = "DNase";
    }

    $("#content" + assembly + assays).show();

    // http://stackoverflow.com/a/8579673
    var aTag = $("a[href='#search']");
    $('html,body').animate({scrollTop: aTag.offset().top},
                           1400,
                           "easeOutQuint");
}

function updateExamples(){
    $(".examples").hide();
    var assembly = activeAssembly()
    $("#examples" + assembly).show();
    $("#lociSearchBox").val(annotationDefaults[assembly]["pos"]);
}

function selectAllNone(setChecked){
    // filter performance note by https://api.jquery.com/hidden-selector/
    $(".selectTissues").filter(":visible")
        .find(":checkbox").prop('checked', setChecked);
}

function select(devPt, checked){
    $(".selectTissues").filter(":visible")
        .find("*[data-selector-name='" + devPt + "']")
        .prop('checked', checked);
}

$(document).ready(function(){
    $('[data-toggle="tooltip"]').tooltip()

    $("#selectAssembly :input").change(function() {
        updateTissues();
        updateExamples();
    });

    // http://stackoverflow.com/a/28603203
    $("#selectAssays :input").change(function() {
        updateTissues();
    });

    $("#selectAll").on('click', function(e){ e.preventDefault(); selectAllNone(true); });
    $("#selectNone").on('click', function(e){ e.preventDefault(); selectAllNone(false); });

    $.each(selectorNames, function(idx, s){
        var klass = ".select" + s;
        $(klass).on('change', function() {
            var val = this.checked;
            select(s, val);
        });
    });

    $("#hg19SNP").click(function(){
        $("#lociSearchBox").val(annotationDefaults["hg19"]["snp"]);
    });
    $("#mm10SNP").click(function(){
        $("#lociSearchBox").val(annotationDefaults["mm10"]["snp"]);
    });
    $("#mm9SNP").click(function(){
        $("#lociSearchBox").val(annotationDefaults["mm9"]["snp"]);
    });

    $('.tablesorter').each(function(i, obj) {
        new Tablesort(obj);
    });

    $( "#searchForm" ).submit(function( event ) {
        event.preventDefault();

        var formData = $("#searchForm").serializeJSON();

        $.ajax({
            type: "POST",
            url: "ucsc",
            data: formData,
            dataType: "json",
            contentType : "application/json",
            success: function(got){
                if("url" in got){
                    return window.open(got["url"], '_blank');
                }
                var w = window.open();
                $(w.document.body).html(got["html"]);
            }
        });
    });

    $("#selectIntersect").click(function() {
        var formData = $("#searchForm").serializeJSON();

        $.ajax({
            type: "POST",
            url: "bedsInRange",
            data: formData,
            dataType: "json",
            contentType : "application/json",
            success: function(got){
                console.log(got);
            }
        });
    });

})

