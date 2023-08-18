

function selectTarget(chembl_id) {
    // Submit form for selected target
    let form = document.getElementById("target_selection");
    let input = document.createElement("input");
    input.type = "hidden"
    input.name = "selected_target";
    input.value = chembl_id;
    form.appendChild(input);
    form.submit();
}

function lipinski() {
    // Submit form to compute Lipinski descriptors
    let form = document.getElementById("lipinski");
    form.submit();
}

function compareModels() {
    // Submit form to compute ML models and compare
    let form = document.getElementById("compareModels");
    form.submit();
}