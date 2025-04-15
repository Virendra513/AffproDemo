
  document.addEventListener("DOMContentLoaded", function() {
    const form = document.getElementById("predict-form");
    const submitButton = document.getElementById("submit-button");
    const otherActionButton = document.getElementById("other-action-button");

    // Handle Submit Button click
    form.addEventListener("submit", function(event) {
        event.preventDefault();  // Prevent form submission initially
        
        // Your form submission logic here
        alert("Form submitted!");

        // Now submit the form
        form.submit();  // If you're ready to submit the form
    });

    // Handle Other Action Button click (for other functionality)
    otherActionButton.addEventListener("click", function() {
        event.preventDefault();  // Prevent form submission initially
        
        // Your form submission logic here
        alert("Form submitted!");

        // Now submit the form
        form.submit();  // If you're ready to submit the form
        // Add your custom logic here for what this button should do (like an API request, show a modal, etc.)
    });
});
