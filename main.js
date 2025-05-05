let imageDates = [];
let index = 0;

function updateImage() {
  document.getElementById("plot").src = imageNames[index];
}

function prevImage() {
  if (index > 0) { index--; updateImage(); }
}

function nextImage() {
  if (index < imageDates.length - 1) { index++; updateImage(); }
}

function loadPage(name) {
document.addEventListener("DOMContentLoaded", () => {
  fetch("file_lists/"+name+".json")
    .then(response => response.json())
    .then(data => {
      imageNames = data;
      index = imageNames.length - 1; // most recent
      updateImage();
    });
});
}
