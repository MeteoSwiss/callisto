let imageDates = [];
let index = 0;

function updateImage() {
  document.getElementById("plot").src = imageNames[index];
}

function prevImage() {
  if (index > 0) { index--; updateImage(); }
}

function nextImage() {
  if (index < imageNames.length - 1) { index++; updateImage(); }
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

window.addEventListener("keydown", event => {
  if(event.keyCode === 39){
  	// next event
    nextImage();
  } else if(event.keyCode === 37){
  	// prev event
    prevImage();
  }
});
